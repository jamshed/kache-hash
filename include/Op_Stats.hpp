
#ifndef KACHE_HASH_OP_STATS_HPP
#define KACHE_HASH_OP_STATS_HPP



#include "utility.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <deque>
#include <mutex>
#include <ostream>
#include <algorithm>
#include <chrono>


namespace kache_hash
{


#ifdef INSTRUMENT
constexpr bool instrument = true;   // Whether the hash table operations are instrumented.
#else
constexpr bool instrument = false;  // Whether the hash table operations are instrumented.
#endif


// Statistics of the hash table operations of a thread, collected iff
// `instrument`. The statistics of each thread are kept in a registry outliving
// the thread, so that the statistics of all the threads, including transient
// ones like resize workers, can be aggregated at quiescent points. An operation
// concludes when the thread's next one begins, or at aggregation, with its
// outcome set by its key getting found or placed. Notes for operations from a
// thread yet to begin any, like a resize worker, are ignored; notes of writer
// locks and of resizes are not.
class Op_Stats
{
public:

    // Types of the operations: updates are insertions and upserts.
    enum Op_Type: uint8_t { update, query, op_type_count };

    // Levels of the table an operation may reach.
    enum Level: uint8_t { primary, secondary, overflow, level_count };

    // Outcomes of the operations: whether the key was found present, placed
    // anew, or found absent.
    enum Outcome: uint8_t { found, placed, absent, outcome_count };

    // Kinds of memory lines touched by the operations: lines of the bucket
    // metadata, of the bucket slots that are first read or first written, and
    // of the overflow table.
    enum Line: uint8_t { meta, slot_read, slot_write, overflow_line, line_count };

private:

    static constexpr std::size_t line_sz = 64;  // Line size for the accounting, fixed to x86's.
    static constexpr std::size_t op_line_cap = 16;  // Maximum number of distinct lines remembered per operation.
    static constexpr std::size_t run_line_cap = 256;    // Maximum number of distinct lines remembered per run.
    static constexpr std::size_t run_len_cap = 256; // Runs at least this long are binned together.

    Op_Type type = update;  // Type of the current operation.
    Level level = primary;  // Deepest level reached by the current operation.
    Outcome outcome = absent;   // Outcome of the current operation so far.
    std::size_t op_line_n = 0;  // Number of distinct lines remembered for the current operation.
    uintptr_t op_line[op_line_cap] = {};    // Lines touched by the current operation.

    uint64_t run_h = 0; // Minimizer hash of the current run: a maximal stretch of operations with the same minimizer hash.
    uint64_t run_len = 0;   // Length of the current run; 0 iff no operation has begun since the last reset.
    std::size_t run_line_n = 0; // Number of distinct lines remembered for the current run.
    uintptr_t run_line[run_line_cap] = {};  // Lines touched by the current run.
    uint64_t run_distinct_c = 0;    // Number of distinct lines touched by the current run.

    time_point_t writer_lock_s = {};    // Beginning of the ongoing writer-lock acquisition.

public:

    uint64_t op_c[op_type_count] = {};  // Number of operations.
    uint64_t concluded_c[op_type_count][level_count][outcome_count] = {};   // Number of operations per deepest level reached and outcome.
    uint64_t line_c[op_type_count][line_count] = {};    // Number of distinct lines touched, summed over the operations.
    uint64_t key_cmp_c[op_type_count] = {}; // Number of key comparisons after checksum matches.
    uint64_t cs_false_c[op_type_count] = {};    // Number of key comparisons after checksum matches with unequal keys.
    uint64_t blind_read_c[op_type_count] = {};  // Number of key reads without checksum matches.
    uint64_t ov_visit_c[op_type_count] = {};    // Number of operations reaching the overflow table.
    uint64_t ov_probe_c[op_type_count] = {};    // Number of slots probed in the overflow table.
    uint64_t ov_lock_c[op_type_count] = {}; // Number of overflow-table slot locks taken.

    uint64_t bucket_lock_c = 0; // Number of compare-and-swap attempts to lock buckets.
    uint64_t bucket_lock_fail_c = 0;    // Number of failed compare-and-swap attempts to lock buckets.
    uint64_t reader_lock_c = 0; // Number of reader-lock acquisitions.
    uint64_t reader_lock_retry_c = 0;   // Number of reader-lock declarations retracted for a writer.
    uint64_t reader_lock_wait_c = 0;    // Number of reader-lock waits for a writer.
    uint64_t writer_lock_c = 0; // Number of writer-lock acquisitions.
    uint64_t writer_lock_wait_ns = 0;   // Time spent acquiring writer locks, in nanoseconds.
    uint64_t resize_c = 0;  // Number of resizes.
    uint64_t resize_ns = 0; // Time spent in resizes, in nanoseconds.
    uint64_t resize_placed_c[level_count] = {}; // Number of keys placed per level of the resized tables.

    uint64_t run_c = 0; // Number of completed runs.
    uint64_t run_len_sum = 0;   // Sum of the completed runs' lengths.
    uint64_t run_len_sq_sum = 0;    // Sum of the completed runs' squared lengths.
    uint64_t run_len_hist[run_len_cap] = {};    // `run_len_hist[i]` is the number of completed runs of length `i + 1`; the last bin is for all longer runs too.
    uint64_t run_len_line_c[run_len_cap] = {};  // `run_len_line_c[i]` is the number of distinct lines touched, summed over the completed runs of `run_len_hist[i]`.
    uint64_t run_line_c[line_count] = {};   // Number of distinct lines touched, summed over the runs.
    uint64_t unremembered_line_c = 0;   // Number of line touches past the remembering capacities, whose lines may be counted again.


private:

    // Returns the registry of the statistics of all the threads.
    static std::deque<Padded<Op_Stats>>& registry();

    // Returns the lock of the registry.
    static std::mutex& registry_lock();

    // Returns whether an operation has begun since the last reset.
    bool in_op() const { return run_len > 0; }

    // Concludes the current operation.
    void end_op() { concluded_c[type][level][outcome]++; }

    // Completes the current run, if any.
    void end_run();

    // Adds the statistics `s` to these, counting its ongoing operation as
    // concluded and its ongoing run as completed.
    void add(const Op_Stats& s);

public:

    // Returns the statistics of the invoking thread.
    static Op_Stats& local();

    // Resets the statistics of all the threads. Not to be invoked concurrently
    // with any operation.
    static void reset();

    // Returns the statistics summed over all the threads. Not to be invoked
    // concurrently with any operation.
    static Op_Stats total();

    // Begins an operation of type `t` for a k-mer with minimizer hash `h`,
    // concluding the previous operation, if any.
    void begin_op(Op_Type t, uint64_t h);

    // Notes that the current operation reached level `l`.
    void reach(Level l);

    // Notes that the current operation touched the memory lines of kind `kind`
    // spanning the `sz` bytes at address `p`. A line counts once per operation
    // and once per run, as the kind of its first touch.
    void touch(Line kind, const void* p, std::size_t sz = 1);

    // Notes a read of the bucket-slot key at address `p` by the current
    // operation, after a checksum match iff `matched`, with the key being equal
    // to the operation's iff `equal`—the key is then found.
    void read_key(const void* p, bool matched, bool equal);

    // Notes that the current operation found its key, accessing the `sz` bytes
    // at address `p` of memory lines of kind `kind`.
    void find_key(Line kind, const void* p, std::size_t sz = 1) { if(in_op()) touch(kind, p, sz), outcome = found; }

    // Notes that the current operation placed its key, writing the `sz` bytes
    // at address `p` of memory lines of kind `kind`.
    void place_key(Line kind, const void* p, std::size_t sz = 1) { if(in_op()) touch(kind, p, sz), outcome = placed; }

    // Notes a probe of the current operation into the overflow-table slot with
    // its key at address `p`.
    void probe_overflow(const void* p) { if(in_op()) touch(overflow_line, p), ov_probe_c[type]++; }

    // Notes an overflow-table slot lock, at address `p`, taken by the current
    // operation.
    void lock_overflow(const void* p) { if(in_op()) touch(overflow_line, p), ov_lock_c[type]++; }

    // Notes a compare-and-swap attempt of the current operation to lock the
    // bucket with its metadata at address `p`.
    void lock_bucket(const void* p) { if(in_op()) touch(meta, p), bucket_lock_c++; }

    // Notes that the latest attempt of the current operation to lock a bucket
    // failed.
    void fail_bucket_lock() { if(in_op()) bucket_lock_fail_c++; }

    // Notes a reader-lock acquisition of the current operation.
    void lock_reader() { if(in_op()) reader_lock_c++; }

    // Notes a wait of the current operation for a writer, to acquire the reader
    // lock.
    void wait_for_writer() { if(in_op()) reader_lock_wait_c++; }

    // Notes that the current operation retracted its reader-lock declaration
    // for a writer.
    void retract_reader() { if(in_op()) reader_lock_retry_c++; }

    // Notes that the invoking thread begins to acquire a writer lock.
    void begin_writer_lock() { writer_lock_s = now(); }

    // Notes that the invoking thread acquired a writer lock.
    void end_writer_lock();

    // Notes a resize by the invoking thread, taking time `d`.
    void resize_table(const std::chrono::nanoseconds d) { resize_c++, resize_ns += d.count(); }

    // Notes that a resize placed a key at level `l` of the new table.
    void reinsert_key(const Level l) { resize_placed_c[l]++; }

    // Writes the statistics to `os` as TSV lines, with keys prefixed by
    // `prefix`.
    void report(std::ostream& os, const std::string& prefix) const;
};


inline std::deque<Padded<Op_Stats>>& Op_Stats::registry()
{
    static std::deque<Padded<Op_Stats>> R;  // A deque keeps the statistics in place as it grows.
    return R;
}


inline std::mutex& Op_Stats::registry_lock()
{
    static std::mutex m;
    return m;
}


inline Op_Stats& Op_Stats::local()
{
    thread_local Op_Stats* const s = []()
    {
        std::lock_guard<std::mutex> guard(registry_lock());
        return &registry().emplace_back().unwrap();
    }();

    return *s;
}


inline void Op_Stats::reset()
{
    std::lock_guard<std::mutex> guard(registry_lock());
    std::for_each(registry().begin(), registry().end(), [](auto& s){ s.unwrap() = Op_Stats(); });
}


inline Op_Stats Op_Stats::total()
{
    std::lock_guard<std::mutex> guard(registry_lock());
    Op_Stats t;
    std::for_each(registry().cbegin(), registry().cend(), [&](const auto& s){ t.add(s.unwrap()); });
    return t;
}


inline void Op_Stats::begin_op(const Op_Type t, const uint64_t h)
{
    if(run_len > 0)
        end_op();

    if(run_len == 0 || h != run_h)
        end_run(), run_h = h;

    run_len++;
    type = t, level = primary, outcome = absent, op_line_n = 0;
    op_c[t]++;
}


inline void Op_Stats::reach(const Level l)
{
    if(in_op() && l > level)
    {
        level = l;
        if(l == overflow)
            ov_visit_c[type]++;
    }
}


inline void Op_Stats::touch(const Line kind, const void* const p, const std::size_t sz)
{
    if(!in_op())
        return;

    const auto addr = reinterpret_cast<uintptr_t>(p);
    for(auto line = addr / line_sz; line <= (addr + sz - 1) / line_sz; ++line)
    {
        if(std::find(op_line, op_line + op_line_n, line) == op_line + op_line_n)
        {
            line_c[type][kind]++;
            if(op_line_n < op_line_cap)
                op_line[op_line_n++] = line;
            else
                unremembered_line_c++;
        }

        if(std::find(run_line, run_line + run_line_n, line) == run_line + run_line_n)
        {
            run_line_c[kind]++, run_distinct_c++;
            if(run_line_n < run_line_cap)
                run_line[run_line_n++] = line;
            else
                unremembered_line_c++;
        }
    }
}


inline void Op_Stats::read_key(const void* const p, const bool matched, const bool equal)
{
    if(!in_op())
        return;

    touch(slot_read, p);
    if(matched)
        key_cmp_c[type]++, cs_false_c[type] += !equal;
    else
        blind_read_c[type]++;

    if(equal)
        outcome = found;
}


inline void Op_Stats::end_run()
{
    if(run_len > 0)
    {
        const auto bin = std::min<uint64_t>(run_len, run_len_cap) - 1;
        run_c++;
        run_len_sum += run_len, run_len_sq_sum += run_len * run_len;
        run_len_hist[bin]++, run_len_line_c[bin] += run_distinct_c;
    }

    run_len = 0, run_line_n = 0, run_distinct_c = 0;
}


inline void Op_Stats::end_writer_lock()
{
    writer_lock_c++;
    writer_lock_wait_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(now() - writer_lock_s).count();
}


inline void Op_Stats::add(const Op_Stats& s)
{
    for(std::size_t t = 0; t < op_type_count; ++t)
    {
        op_c[t] += s.op_c[t];
        for(std::size_t l = 0; l < level_count; ++l)
            for(std::size_t o = 0; o < outcome_count; ++o)
                concluded_c[t][l][o] += s.concluded_c[t][l][o];

        for(std::size_t kind = 0; kind < line_count; ++kind)
            line_c[t][kind] += s.line_c[t][kind];

        key_cmp_c[t] += s.key_cmp_c[t], cs_false_c[t] += s.cs_false_c[t], blind_read_c[t] += s.blind_read_c[t];
        ov_visit_c[t] += s.ov_visit_c[t], ov_probe_c[t] += s.ov_probe_c[t], ov_lock_c[t] += s.ov_lock_c[t];
    }

    bucket_lock_c += s.bucket_lock_c, bucket_lock_fail_c += s.bucket_lock_fail_c;
    reader_lock_c += s.reader_lock_c, reader_lock_retry_c += s.reader_lock_retry_c, reader_lock_wait_c += s.reader_lock_wait_c;
    writer_lock_c += s.writer_lock_c, writer_lock_wait_ns += s.writer_lock_wait_ns;
    resize_c += s.resize_c, resize_ns += s.resize_ns;
    for(std::size_t l = 0; l < level_count; ++l)
        resize_placed_c[l] += s.resize_placed_c[l];

    run_c += s.run_c, run_len_sum += s.run_len_sum, run_len_sq_sum += s.run_len_sq_sum;
    for(std::size_t i = 0; i < run_len_cap; ++i)
        run_len_hist[i] += s.run_len_hist[i], run_len_line_c[i] += s.run_len_line_c[i];

    for(std::size_t kind = 0; kind < line_count; ++kind)
        run_line_c[kind] += s.run_line_c[kind];

    unremembered_line_c += s.unremembered_line_c;

    if(s.run_len > 0)   // The ongoing operation and run of `s` count as concluded and completed; their lines are already counted.
    {
        const auto bin = std::min<uint64_t>(s.run_len, run_len_cap) - 1;
        concluded_c[s.type][s.level][s.outcome]++;
        run_c++;
        run_len_sum += s.run_len, run_len_sq_sum += s.run_len * s.run_len;
        run_len_hist[bin]++, run_len_line_c[bin] += s.run_distinct_c;
    }
}


inline void Op_Stats::report(std::ostream& os, const std::string& prefix) const
{
    constexpr const char* type_name[op_type_count] = {"update", "query"};
    constexpr const char* level_name[level_count] = {"primary", "secondary", "overflow"};
    constexpr const char* outcome_name[outcome_count] = {"found", "placed", "absent"};
    constexpr const char* line_name[line_count] = {"meta_lines", "slot_read_lines", "slot_write_lines", "overflow_lines"};

    const auto put = [&](const std::string& key, const uint64_t val){ os << prefix << key << "\t" << val << "\n"; };

    for(std::size_t t = 0; t < op_type_count; ++t)
    {
        const std::string ty(type_name[t]);
        put(ty + ".ops", op_c[t]);
        for(std::size_t l = 0; l < level_count; ++l)
            for(std::size_t o = 0; o < outcome_count; ++o)
                put(ty + "." + level_name[l] + "." + outcome_name[o], concluded_c[t][l][o]);

        for(std::size_t kind = 0; kind < line_count; ++kind)
            put(ty + "." + line_name[kind], line_c[t][kind]);

        put(ty + ".key_compares", key_cmp_c[t]);
        put(ty + ".checksum_false_matches", cs_false_c[t]);
        put(ty + ".blind_reads", blind_read_c[t]);
        put(ty + ".overflow_visits", ov_visit_c[t]);
        put(ty + ".overflow_probes", ov_probe_c[t]);
        put(ty + ".overflow_locks", ov_lock_c[t]);
    }

    put("bucket_lock.tries", bucket_lock_c);
    put("bucket_lock.fails", bucket_lock_fail_c);
    put("reader_lock.acquisitions", reader_lock_c);
    put("reader_lock.retries", reader_lock_retry_c);
    put("reader_lock.waits", reader_lock_wait_c);
    put("writer_lock.acquisitions", writer_lock_c);
    put("writer_lock.wait_ns", writer_lock_wait_ns);
    put("resize.count", resize_c);
    put("resize.ns", resize_ns);
    for(std::size_t l = 0; l < level_count; ++l)
        put(std::string("resize.placed.") + level_name[l], resize_placed_c[l]);

    put("lines.unremembered", unremembered_line_c);
    put("runs.count", run_c);
    put("runs.len_sum", run_len_sum);
    put("runs.len_sq_sum", run_len_sq_sum);
    for(std::size_t kind = 0; kind < line_count; ++kind)
        put(std::string("runs.") + line_name[kind], run_line_c[kind]);

    for(std::size_t i = 0; i < run_len_cap; ++i)
        if(run_len_hist[i] > 0)
        {
            const auto len = std::to_string(i + 1) + (i + 1 == run_len_cap ? "+" : "");
            put("runs.len." + len, run_len_hist[i]);
            put("runs.lines." + len, run_len_line_c[i]);
        }
}


}



#endif
