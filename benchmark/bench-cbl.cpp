
#include "bench.hpp"

#include <cstddef>
#include <cstdlib>
#include <iostream>


using namespace kache_hash;


// C interface to the CBL set, from the crate in `cbl/`.
extern "C"
{
    struct CBL_Set; // Opaque CBL set.

    CBL_Set* cbl_new();
    void cbl_free(CBL_Set* s);
    void cbl_insert_seq(CBL_Set* s, const char* seq, std::size_t len);
    std::size_t cbl_count_seq(CBL_Set* s, const char* seq, std::size_t len);
    std::size_t cbl_count(const CBL_Set* s);
    double cbl_prefix_load(const CBL_Set* s);
    std::size_t cbl_k();
    std::size_t cbl_prefix_bits();
}


// Adapter of the CBL set—sequential, and with no notion of capacity.
class CBL
{
private:

    CBL_Set* const S;   // The set.

public:

    static constexpr auto name = "cbl"; // Name of the structure.
    static constexpr bool is_set = true;    // Whether a set.
    static constexpr bool concurrent = false;   // Whether safe for multiple user threads.

    // Prints the parameters of the set.
    static void report_params() { print("cbl.prefix_bits", cbl_prefix_bits()); }

    // Constructs an empty set; CBL can not be sized in advance.
    CBL(std::size_t, std::size_t): S(cbl_new()) {}

    CBL(const CBL&) = delete;
    CBL& operator=(const CBL&) = delete;

    ~CBL() { cbl_free(S); }

    // Adds the `n` k-mers of the sequence `s`.
    void insert(std::size_t, const char* s, std::size_t n) { cbl_insert_seq(S, s, n + k - 1); }

    // Returns the number of the `n` k-mers of the sequence `s` present in the
    // set.
    std::size_t find(std::size_t, const char* s, std::size_t n) const { return cbl_count_seq(S, s, n + k - 1); }

    // Returns the number of k-mers in the set.
    std::size_t size() const { return cbl_count(S); }

    // Prints the statistics of the set.
    void report() const { print("cbl.prefix_load", cbl_prefix_load(S)); }
};


int main(int argc, char* argv[])
{
    const auto cfg = parse_config(argc, argv);

    if(cbl_k() != k)
    {
        std::cerr << "CBL is built for k = " << cbl_k() << ", and the driver for k = " << k << ".\n";
        std::exit(EXIT_FAILURE);
    }

    if(!cfg.set)
        reject("CBL is a set, and hence can not count");

    if(!cfg.seq)
        reject("CBL has no concurrent variant");

    bench<CBL>(cfg);

    return 0;
}
