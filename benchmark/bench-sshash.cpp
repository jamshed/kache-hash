
#include "bench.hpp"

#include "include/dictionary_types.hpp"
#include "include/streaming_query.hpp"
#include "src/dictionary.cpp"   // SSHash's templates are defined in its sources, which its tools include likewise.
#include "src/query.cpp"

#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>


using namespace kache_hash;


// Adapter of the SSHash dictionary: a static set of k-mers, loaded from its
// index file and looked up through a streaming query per user thread.
class SSHash
{
private:

    typedef sshash::dictionary_type dict_t;
    typedef sshash::streaming_query<dict_t> query_t;

    dict_t D;   // The dictionary.
    mutable std::vector<Padded<query_t>> Q; // `Q[t]` is the streaming query of the `t`'th user thread.

public:

    static constexpr auto name = "sshash";  // Name of the structure.
    static constexpr bool is_set = true;    // Whether a set.
    static constexpr bool concurrent = true;    // Whether safe for multiple user threads.
    static constexpr bool is_static = true; // Whether loaded from an index file rather than built.

    // Prints the parameters of the dictionary.
    static void report_params() {}

    // Loads the dictionary from the index file at path `path`, for `thread_c`
    // user threads.
    SSHash(const std::string& path, std::size_t thread_c);

    // The streaming queries point into the dictionary.
    SSHash(const SSHash&) = delete;
    SSHash& operator=(const SSHash&) = delete;

    // Returns the number of the `n` k-mers of the sequence `s` present in the
    // dictionary, for the `t_id`'th user thread.
    std::size_t find(std::size_t t_id, const char* s, std::size_t n) const;

    // Returns the number of k-mers in the dictionary.
    std::size_t size() const { return D.num_kmers(); }

    // Prints the statistics of the dictionary.
    void report() const;
};


inline SSHash::SSHash(const std::string& path, const std::size_t thread_c)
{
    essentials::load(D, path.c_str());
    if(D.k() != k)
    {
        std::cerr << "The SSHash index is built for k = " << D.k() << ", and the driver for k = " << k << ".\n";
        std::exit(EXIT_FAILURE);
    }

    Q.reserve(thread_c);
    for(std::size_t t = 0; t < thread_c; ++t)
        Q.emplace_back(query_t(&D));
}


inline std::size_t SSHash::find(const std::size_t t_id, const char* const s, const std::size_t n) const
{
    auto& q = Q[t_id].unwrap();
    q.reset();  // The run does not continue the last one.

    std::size_t h = 0;
    for(std::size_t i = 0; i < n; ++i)
        h += (q.lookup(s + i).kmer_id != sshash::constants::invalid_uint64);

    return h;
}


inline void SSHash::report() const
{
    print("sshash.m", D.m());
    print("sshash.num_strings", D.num_strings());
    print("sshash.bits_per_kmer", static_cast<double>(D.num_bits()) / D.num_kmers());
}


int main(int argc, char* argv[])
{
    const auto cfg = parse_config(argc, argv);

    if(!cfg.set)
        reject("SSHash is a set, and hence can not count");

    if(cfg.seq)
        reject("SSHash has no sequential variant");

    bench<SSHash>(cfg);

    return 0;
}
