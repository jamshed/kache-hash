#include "Streaming_Kmer_Hash_Table.hpp"
#include "Parser.hpp"

#include <iostream>


void hash_set_sequential()
{
    constexpr uint16_t k = 31;  // k-mer size.
    constexpr uint16_t l = 19;  // Minimizer size.

    // Hash set with an initial capacity for 100M k-mers, that will use 1 worker thread during resizes.
    kache_hash::Streaming_Kmer_Hash_Table<k, false> S(100'000'000, 1);

    kache_hash::Parser input("data/ecoli.fasta");   // Parser to read through the input file.
    const auto token = S.register_user();   // The user thread needs to be registered to use the set and get a token to use with all the operations.

    while(input.read_next_seq())    // Read the input sequences from the file
    {
        auto const seq = input.seq();   // Each sequence.
        auto const len = input.seq_len();   // Length of the sequence.
        if(len < k)
            continue;

        kache_hash::Kmer_Window<k, l> w;    // Window to stream through the k-mers of the sequence.
        w.init(seq);    // Initialize the k-mer window at the first k-mer of `seq`.
        std::size_t kmer_idx = 0;   // Index of the current k-mer at sequence `seq`.
        while(true)
        {
            S.insert(w, token); // The user thread having the token `token` inserts the k-mer defined by the window `w` into the set.
            if(kmer_idx + k == len)
                break;

            w.advance(seq[kmer_idx + k]);   // Move the window to the next k-mer.
            kmer_idx++;
        }
    }

    std::cout << "#Unique k-mers: " << S.size() << "\n";
}


int main()
{
    hash_set_sequential();

    return 0;
}
