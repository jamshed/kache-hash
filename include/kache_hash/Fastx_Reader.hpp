
#ifndef KACHE_HASH_FASTX_READER_HPP
#define KACHE_HASH_FASTX_READER_HPP



#include <cstddef>
#include <functional>
#include <memory>
#include <string>


// Marks the library's classes to export from a shared library, which hides
// its other symbols.
#ifndef KACHE_COUNTER_EXPORT
#define KACHE_COUNTER_EXPORT __attribute__((visibility("default")))
#endif


namespace kache_hash
{

// Parallel reader of the sequences of a FASTA or FASTQ file; interface designed
// so that it can directly feed Jellyfish. The first character sets the format:
// '>' for FASTA, '@' for FASTQ. A FASTA sequence joins the lines up to the next
// header; a FASTQ record has four lines, and its quality line is as long as its
// sequence. Lines end with "\n" or "\r\n". A regular file is memory-mapped and
// split among the workers; any other input, such as a pipe, is read by one more
// thread, which hands them chunks of whole records.
class KACHE_COUNTER_EXPORT Fastx_Reader
{
private:

    struct Impl;
    std::unique_ptr<Impl> impl;

public:

    // Callback of the worker `t` with a batch of the `c` sequences
    // `seqs[0..c)`, of lengths `lens[0..c)`.
    typedef std::function<void(std::size_t t, const char* const* seqs, const std::size_t* lens, std::size_t c)> batch_f;

    // Callback of the worker `t` after its last batch.
    typedef std::function<void(std::size_t t)> done_f;

    // Bytes of input per unit of work.
    static constexpr std::size_t chunk_bytes_default = std::size_t(16) << 20;

    // Constructs a reader of the file at `path`, which hands out its input in
    // units of about `chunk_bytes` bytes.
    explicit Fastx_Reader(const std::string& path, std::size_t chunk_bytes = chunk_bytes_default);

    ~Fastx_Reader();

    // Calls `f` on `thread_c` workers with the input's non-empty sequences,
    // each once, in an unspecified order, in batches of at most 1024 sequences
    // and about 1 Mbp; then `done(t)` on each worker `t`, if `done` is given.
    // Throws `std::runtime_error` on malformed input, or what the callbacks
    // throw. The input is read once: a reader reads only once.
    void read(std::size_t thread_c, const batch_f& f, const done_f& done = nullptr);
};

}



#endif
