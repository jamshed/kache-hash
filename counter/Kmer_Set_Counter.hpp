
#ifndef KACHE_HASH_KMER_SET_COUNTER_HPP
#define KACHE_HASH_KMER_SET_COUNTER_HPP



#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>


namespace kache_hash
{

// Counter of the k-mers of reads that belong to a set fixed beforehand, as
// PanGenie counts them with Jellyfish:  K-mers are canonical: a k-mer and its
// reverse complement are treated as one. Inputs are FASTA or FASTQ files, and
// non-ACGT characters split their sequences. Counts are 32-bit.
// The counter's operations must not overlap one another, but for `histogram`s
// and `count_of`s.
class Kmer_Set_Counter
{
private:

    struct Impl;
    std::unique_ptr<Impl> impl;

public:

    // Time taken by the phases, in seconds.
    struct Timings
    {
        double alloc = 0;   // Mapping the table.
        double prime = 0;
        double count = 0;
        double histogram = 0;
    };

    // Returns whether k-mers of length `k` are supported: only 31-mers are, for
    // now.
    static bool supports(uint16_t k);

    // Constructs a counter of `k`-mers with room for `n` k-mers at load factor
    // `lf`, beyond which it grows, with its table's bucket count rounded up to
    // a power of 2, and with `thread_c` threads in `[1, 128]` for its
    // operations.
    Kmer_Set_Counter(uint16_t k, std::size_t n, std::size_t thread_c, double lf = 0.9);

    ~Kmer_Set_Counter();

    // Adds the k-mers of the sequences of the FASTA or FASTQ file at `path` to
    // the set, with count 0; the ones present keep their counts.
    void prime(const std::string& path);

    // Increments the count of each k-mer of the reads of the FASTA or FASTQ
    // file at `path` that is in the set, at each of its occurrences. Returns
    // the number of occurrences counted.
    std::size_t count(const std::string& path);

    // Returns the histogram of the counts: the number of k-mers with count `c`
    // at index `c`, for `c` in `[1, max]`, with 0 at index 0; zero and larger
    // counts are left out.
    std::vector<std::size_t> histogram(std::size_t max);

    // Returns the count of the k-mer of the `k` characters at `kmer`: 0 if it
    // is not in the set, as for a k-mer with a non-ACGT character.
    uint64_t count_of(const char* kmer) const;

    // Returns the number of k-mers in the set.
    std::size_t size() const;

    // Returns the time taken by the phases so far.
    const Timings& timings() const;
};


}



#endif
