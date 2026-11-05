
#include "Fastx_Reader.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <condition_variable>
#include <cstring>
#include <deque>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>


namespace kache_hash
{

namespace
{

constexpr std::size_t batch_seqs = 1024;    // Most sequences in a batch.
constexpr std::size_t batch_bases = std::size_t(1) << 20;   // Bases at which a batch is handed over.
constexpr std::size_t cut_window = std::size_t(1) << 20;    // Tail of a stream's buffer where its cut is looked for first.

enum class Format { fasta, fastq };


// Returns the format of the input `path`, whose first bytes are `p[0..n)`,
// `n > 0`.
Format format_of(const char* const p, const std::size_t n, const std::string& path)
{
    if(n >= 2 && static_cast<unsigned char>(p[0]) == 0x1f && static_cast<unsigned char>(p[1]) == 0x8b)
        throw std::runtime_error(path + " is gzip-compressed; decompress it, e.g. into a pipe");
    if(p[0] == '>')
        return Format::fasta;
    if(p[0] == '@')
        return Format::fastq;

    throw std::runtime_error(path + " is neither FASTA nor FASTQ: it starts with neither '>' nor '@'");
}


// Returns the offset of the first record of the `format` text `p[0..n)` that
// starts at or after `off`, or `n` if none does. A FASTA record starts at a
// line beginning with '>'. A FASTQ record starts at a line beginning with '@'
// whose line two after begins with '+': a quality line beginning with '@' is
// followed by a header, and then a sequence.
std::size_t record_start(const char* const p, const std::size_t n, std::size_t off, const Format format)
{
    if(off == 0)
        return 0;

    while(off < n)
    {
        if(format == Format::fasta)
        {
            const char* const r = static_cast<const char*>(std::memchr(p + off, '>', n - off));
            if(!r)
                return n;
            if(r[-1] == '\n')
                return r - p;

            off = r - p + 1;
            continue;
        }

        const char* const nl = static_cast<const char*>(std::memchr(p + off - 1, '\n', n - off + 1));
        if(!nl || nl + 1 == p + n)
            return n;

        const std::size_t s = nl - p + 1;   // The first line start at or after `off`.
        if(p[s] == '@')
        {
            const char* const l1 = static_cast<const char*>(std::memchr(p + s, '\n', n - s));
            const char* const l2 = l1 ? static_cast<const char*>(std::memchr(l1 + 1, '\n', p + n - l1 - 1)) : nullptr;
            if(l2 && l2 + 1 < p + n && l2[1] == '+')
                return s;
        }

        off = s + 1;
    }

    return n;
}


// A line of a text: the end of its content, before its "\n" or "\r\n", and the
// start of the next line.
struct Line
{
    const char* end;
    const char* next;
};


// Returns the line at `s` of the text ending at `e`.
Line line_at(const char* const s, const char* const e)
{
    const char* const nl = static_cast<const char*>(std::memchr(s, '\n', e - s));
    const char* const end = (nl ? nl : e);
    return {(end > s && end[-1] == '\r') ? end - 1 : end, nl ? nl + 1 : e};
}


// A worker's batch of sequences under construction. A sequence of several lines
// is joined in the batch's buffer.
class Batch
{
private:

    const Fastx_Reader::batch_f& f;
    const std::size_t t;    // The worker.

    std::vector<const char*> seqs;
    std::vector<std::size_t> lens;
    std::vector<std::pair<std::size_t, std::size_t>> joined;    // The joined sequences' indices, and offsets in `buf`.
    std::string buf;
    std::size_t bases = 0;

    // Notes the addition of a sequence of length `len`, and hands over the batch once full.
    void added(const std::size_t len)
    {
        bases += len;
        if(seqs.size() == batch_seqs || bases >= batch_bases)
            flush();
    }

public:

    Batch(const Fastx_Reader::batch_f& f, const std::size_t t):
          f(f)
        , t(t)
    {
        seqs.reserve(batch_seqs), lens.reserve(batch_seqs);
    }

    // Adds the sequence `s[0..len)`, if not empty.
    void add(const char* const s, const std::size_t len)
    {
        if(len > 0)
            seqs.push_back(s), lens.push_back(len), added(len);
    }

    // Returns the offset in the buffer at which a sequence to join would start.
    std::size_t join_off() const { return buf.size(); }

    // Appends the line `s[0..len)` to the sequence being joined.
    void append(const char* const s, const std::size_t len) { buf.append(s, len); }

    // Adds the sequence joined from offset `off` of the buffer, if not empty.
    void add_joined(const std::size_t off)
    {
        const std::size_t len = buf.size() - off;
        if(len > 0)
            joined.emplace_back(seqs.size(), off), seqs.push_back(nullptr), lens.push_back(len), added(len);
    }

    // Hands the batch to the callback, if it is not empty, and clears it.
    void flush()
    {
        if(seqs.empty())
            return;

        for(const auto& [i, off] : joined)
            seqs[i] = buf.data() + off;
        f(t, seqs.data(), lens.data(), seqs.size());
        seqs.clear(), lens.clear(), joined.clear(), buf.clear(), bases = 0;
    }
};


// Adds the sequences of the `format` records of `p[b..e)` to `batch`: `b` is
// the start of a record, and `e` that of another, or the end of the input. `p`
// is at byte `p_off` of the input.
void parse(const char* const p, const std::size_t b, const std::size_t e, const std::size_t p_off, const Format format, Batch& batch)
{
    const char* r = p + b;
    const char* const end = p + e;
    const auto malformed = [&](const char* const what)
    {
        throw std::runtime_error(std::string("malformed ") + (format == Format::fasta ? "FASTA" : "FASTQ") +
                                 " record at byte " + std::to_string(p_off + (r - p)) + ": " + what);
    };

    if(format == Format::fasta)
        while(r < end)
        {
            if(*r != '>')
                malformed("expected a header");

            const char* const s = line_at(r, end).next; // The sequence's first line.
            if(s == end || *s == '>')
            {
                r = s;
                continue;
            }

            const auto L = line_at(s, end);
            if(L.next == end || *L.next == '>')
            {
                batch.add(s, L.end - s), r = L.next;
                continue;
            }

            const auto off = batch.join_off();
            for(r = s; r < end && *r != '>'; )
            {
                const auto M = line_at(r, end);
                batch.append(r, M.end - r);
                r = M.next;
            }

            batch.add_joined(off);
        }
    else
        while(r < end)
        {
            if(*r != '@')
                malformed("expected a header");

            const char* const s = line_at(r, end).next; // The sequence.
            const auto S = line_at(s, end);
            if(S.next == end || *S.next != '+')
                malformed("expected a '+' line");

            const char* const q = line_at(S.next, end).next;    // The quality.
            const auto Q = line_at(q, end);
            if(Q.end - q != S.end - s)
                malformed("its quality and sequence lengths differ");

            batch.add(s, S.end - s);
            r = Q.next;
        }
}


// Reads from the stream `fd` into `buf[filled..want)`, growing `buf` to `want`,
// until it is full or the stream ends. Returns whether it is full.
bool fill(const int fd, std::vector<char>& buf, std::size_t& filled, const std::size_t want)
{
    if(buf.size() < want)
        buf.resize(want);

    while(filled < want)
    {
        const auto r = ::read(fd, buf.data() + filled, want - filled);
        if(r < 0)
        {
            if(errno == EINTR)
                continue;
            throw std::runtime_error(std::string("cannot read the input: ") + std::strerror(errno));
        }

        if(r == 0)
            return false;
        filled += r;
    }

    return true;
}


// A unit of whole records read from a stream: its buffer, its size, and its
// offset in the stream.
struct Unit
{
    std::vector<char> buf;
    std::size_t size = 0;
    std::size_t off = 0;
};


// Bounded queue of units, from the stream's reader to the workers, which also
// keeps the buffers of the units read.
class Unit_Queue
{
private:

    std::mutex m;
    std::condition_variable not_full, not_empty;
    std::deque<Unit> Q;
    std::vector<std::vector<char>> bufs;    // Buffers for reuse.
    std::size_t cap = 1;
    bool closed = false;
    bool aborted = false;

public:

    // Sets the capacity of the queue to `cap` units.
    void set_cap(const std::size_t cap) { std::lock_guard lock(m); this->cap = cap; }

    // Adds the unit `u`, once there is room. Returns whether the queue runs
    // still.
    bool put(Unit&& u)
    {
        std::unique_lock lock(m);
        not_full.wait(lock, [&](){ return Q.size() < cap || aborted; });
        if(aborted)
            return false;

        Q.push_back(std::move(u));
        not_empty.notify_one();
        return true;
    }

    // Takes the next unit into `u`, once there is one. Returns `false` once the
    // queue is closed and empty, or aborted.
    bool get(Unit& u)
    {
        std::unique_lock lock(m);
        not_empty.wait(lock, [&](){ return !Q.empty() || closed || aborted; });
        if(aborted || Q.empty())
            return false;

        u = std::move(Q.front());
        Q.pop_front();
        not_full.notify_one();
        return true;
    }

    // Marks the end of the units.
    void close() { std::lock_guard lock(m); closed = true; not_empty.notify_all(); }

    // Stops the queue: its waits end, and no more units go through.
    void abort() { std::lock_guard lock(m); aborted = true; not_full.notify_all(), not_empty.notify_all(); }

    // Returns a buffer for reuse, or an empty one.
    std::vector<char> take_buf()
    {
        std::lock_guard lock(m);
        if(bufs.empty())
            return {};

        auto b = std::move(bufs.back());
        bufs.pop_back();
        return b;
    }

    // Keeps the buffer `b` for reuse.
    void give_buf(std::vector<char>&& b) { std::lock_guard lock(m); bufs.push_back(std::move(b)); }
};


}


struct Fastx_Reader::Impl
{
    const std::string path;
    const std::size_t chunk_bytes;
    int fd = -1;    // The input, if it is not a regular file.
    const char* p = nullptr;    // The input, if it is a regular non-empty file, mapped.
    std::size_t n = 0;  // Its size.
    bool used = false;  // Whether the input is read.

    std::mutex err_m;
    std::exception_ptr err; // The first exception of the threads.
    std::atomic<bool> stop{false};  // Whether the threads should stop.
    Unit_Queue q;   // The units of a stream.

    Impl(const std::string& path, std::size_t chunk_bytes);

    ~Impl();

    // Records the current exception, unless an earlier one is, and stops the
    // threads.
    void fail();

    // Runs `work(t)` and then `done(t)`, if given, on workers `t` in `[0,
    // thread_c)`, and `extra()` on one more thread, if given. Rethrows the
    // first exception of any.
    void run(std::size_t thread_c, const std::function<void(std::size_t)>& work, const done_f& done, const std::function<void()>& extra = nullptr);

    // Reads the mapped input, as `read` does.
    void read_mapped(std::size_t thread_c, const batch_f& f, const done_f& done);

    // Reads the stream input, as `read` does.
    void read_stream(std::size_t thread_c, const batch_f& f, const done_f& done);
};


Fastx_Reader::Impl::Impl(const std::string& path, const std::size_t chunk_bytes):
      path(path)
    , chunk_bytes(std::max<std::size_t>(chunk_bytes, 1))
{
    fd = ::open(path.c_str(), O_RDONLY);
    if(fd < 0)
        throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));

    struct stat st;
    if(fstat(fd, &st) != 0 || !S_ISREG(st.st_mode))
        return;

    n = st.st_size;
    if(n > 0)
    {
        void* const m = mmap(nullptr, n, PROT_READ, MAP_PRIVATE, fd, 0);
        if(m == MAP_FAILED)
        {
            const auto e = errno;
            ::close(fd);
            throw std::runtime_error("cannot map " + path + ": " + std::strerror(e));
        }

        p = static_cast<const char*>(m);
    }

    ::close(fd), fd = -1;
}


Fastx_Reader::Impl::~Impl()
{
    if(p)
        munmap(const_cast<char*>(p), n);
    if(fd >= 0)
        ::close(fd);
}


void Fastx_Reader::Impl::fail()
{
    {
        std::lock_guard lock(err_m);
        if(!err)
            err = std::current_exception();
    }

    stop = true;
    q.abort();
}


void Fastx_Reader::Impl::run(const std::size_t thread_c, const std::function<void(std::size_t)>& work, const done_f& done, const std::function<void()>& extra)
{
    std::vector<std::thread> W;
    for(std::size_t t = 0; t < thread_c; ++t)
        W.emplace_back([&, t]()
        {
            try { work(t); } catch(...) { fail(); }

            if(done)    // Even after a failure: a table may wait for its threads to be done, e.g. Jellyfish's to grow.
                try { done(t); } catch(...) { fail(); }
        });

    if(extra)
        W.emplace_back([&](){ try { extra(); } catch(...) { fail(); } });

    for(auto& w : W)
        w.join();

    if(err)
        std::rethrow_exception(err);
}


void Fastx_Reader::Impl::read_mapped(const std::size_t thread_c, const batch_f& f, const done_f& done)
{
    const auto format = (n > 0 ? format_of(p, n, path) : Format::fasta);
    const std::size_t page = sysconf(_SC_PAGESIZE);
    std::atomic<std::size_t> next{0};   // The next chunk of the input.

    run(thread_c, [&](const std::size_t t)
    {
        Batch batch(f, t);
        while(!stop)
        {
            const std::size_t a = next++ * chunk_bytes;
            if(a >= n)
                break;

            // The records that start in the chunk.
            const auto b = record_start(p, n, a, format), e = record_start(p, n, std::min(a + chunk_bytes, n), format);
            if(b >= e)
                continue;

            parse(p, b, e, 0, format, batch);
            batch.flush();

            // The pages read stay in the page cache, but leave the process's resident set.
            const auto u = (b + page - 1) / page * page, v = e / page * page;
            if(u < v)
                madvise(const_cast<char*>(p) + u, v - u, MADV_DONTNEED);
        }
    }, done);
}


void Fastx_Reader::Impl::read_stream(const std::size_t thread_c, const batch_f& f, const done_f& done)
{
    // The first chunk sets the format.
    Unit cur;
    bool more = fill(fd, cur.buf, cur.size, chunk_bytes);
    const auto format = (cur.size > 0 ? format_of(cur.buf.data(), cur.size, path) : Format::fasta);
    q.set_cap(2 * thread_c);

    // Hands out the stream in units cut at records' starts: in the tail of the buffer, or else past its first record.
    const auto produce = [&]()
    {
        while(!stop)
        {
            if(!more)
            {
                if(cur.size > 0)
                    q.put(std::move(cur));
                break;
            }

            const char* const b = cur.buf.data();
            auto cut = record_start(b, cur.size, cur.size - std::min(cur.size / 2, cut_window), format);
            if(cut == cur.size)
                cut = record_start(b, cur.size, 1, format);
            if(cut == cur.size)
            {
                more = fill(fd, cur.buf, cur.size, 2 * cur.size);
                continue;
            }

            Unit next;
            next.buf = q.take_buf(), next.size = cur.size - cut, next.off = cur.off + cut;
            if(next.buf.size() < next.size)
                next.buf.resize(next.size);
            std::memcpy(next.buf.data(), b + cut, next.size);

            cur.size = cut;
            if(!q.put(std::move(cur)))
                break;

            cur = std::move(next);
            more = fill(fd, cur.buf, cur.size, std::max(chunk_bytes, 2 * cur.size));
        }

        q.close();
    };

    run(thread_c, [&](const std::size_t t)
    {
        Batch batch(f, t);
        for(Unit u; q.get(u); )
        {
            parse(u.buf.data(), 0, u.size, u.off, format, batch);
            batch.flush();
            q.give_buf(std::move(u.buf));
        }
    }, done, produce);
}


Fastx_Reader::Fastx_Reader(const std::string& path, const std::size_t chunk_bytes):
    impl(std::make_unique<Impl>(path, chunk_bytes))
{}


Fastx_Reader::~Fastx_Reader() = default;


void Fastx_Reader::read(const std::size_t thread_c, const batch_f& f, const done_f& done)
{
    if(impl->used)
        throw std::logic_error("a Fastx_Reader reads its input only once");
    impl->used = true;

    const auto t_c = std::max<std::size_t>(thread_c, 1);
    impl->fd < 0 ? impl->read_mapped(t_c, f, done) : impl->read_stream(t_c, f, done);
}


}
