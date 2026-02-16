#include "Streaming_Kmer_Hash_Table.hpp"
#include "Parser.hpp"

#include <algorithm>
#include <iostream>
#include <thread>
#include <vector>


void hash_set_sequential()
{
    constexpr uint16_t k = 31;  // k-mer size.

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

        kache_hash::Kmer_Window<k> w;    // Window to stream through the k-mers of the sequence.
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


void hash_set_concurrent()
{
    constexpr uint16_t k = 31;  // k-mer size.
    const uint64_t thread_count = 16;   // Number of user threads.

    typedef kache_hash::Streaming_Kmer_Hash_Table<k, true> ht_t;

    // Hash set with an initial capacity for 100M k-mers, that will use `thread_count` worker threads during resizes.
    ht_t S(100'000'000, thread_count);

    std::vector<typename ht_t::Token> token(thread_count);  // Each separate user thread needs a token to use with all hash table operations.
    for(uint16_t w = 0; w < token.size(); ++w)
        token[w] = S.register_user();   // Each user thread needs to be registered to use the set.

    kache_hash::Parser input("data/ecoli.fasta");   // Parser to read through the input file.

    while(input.read_next_seq())    // Read the input sequences from the file.
    {
        auto const seq = input.seq();   // Each sequence.
        auto const len = input.seq_len();   // Length of the sequence.
        if(len < k)
            continue;

        std::vector<std::thread> W; // Worker threads.
        const auto kmer_per_w = (len - k + 1) / thread_count;   // Number of k-mers to be processed from `seq` by each worker thread.
        for(uint64_t t = 0; t < thread_count; ++t)
            // Add each thread
            W.emplace_back([&](const uint64_t t_id)
            {
                const auto s = seq + t_id * kmer_per_w; // Starting point of the sequence to be processed by the `t_id`'th thread.
                const auto s_len = (t_id < thread_count - 1 ? kmer_per_w + k - 1 : len - t_id * kmer_per_w);    // Length of the sequence to be processed by this thread.
                if(s_len < k)
                    return;

                kache_hash::Kmer_Window<k> w;   // Window to stream through the k-mers of the sequence.
                w.init(s);  // Initialize the k-mer window at the first k-mer of the sequence assigned to this thread.
                std::size_t kmer_idx = 0;   // Index of the current k-mer.
                while(true)
                {
                    S.insert(w, token[t_id]);   // The user thread having the `t_id`'th token inserts the k-mer defined by the window `w` into the set.
                    if(kmer_idx + k == s_len)
                        break;

                    w.advance(s[kmer_idx + k]); // Move the window to the next k-mer.
                    kmer_idx++;
                }
            }
            , t);

        std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });    // Join the worker threads.
    }

    std::cout << "#Unique k-mers: " << S.size() << "\n";
}


void hash_map_sequential()
{
    constexpr uint16_t k = 31;  // k-mer size.

    // Hash map for (k-mer, count) key-value pairs with an initial capacity for 100M k-mers, that will use 1 worker thread during resizes.
    kache_hash::Streaming_Kmer_Hash_Table<k, false, uint16_t> M(100'000'000, 1);

    kache_hash::Parser input("data/ecoli.fasta");   // Parser to read through the input file.
    const auto token = M.register_user();   // The user thread needs to be registered to use the map and get a token to use with all the operations.

    while(input.read_next_seq())    // Read the input sequences from the file
    {
        auto const seq = input.seq();   // Each sequence.
        auto const len = input.seq_len();   // Length of the sequence.
        if(len < k)
            continue;

        kache_hash::Kmer_Window<k> w;    // Window to stream through the k-mers of the sequence.
        w.init(seq);    // Initialize the k-mer window at the first k-mer of `seq`.
        std::size_t kmer_idx = 0;   // Index of the current k-mer at sequence `seq`.
        while(true)
        {
            // The user thread having the token `token` either inserts the k-mer defined by the window `w` into the map with a count of 1 if it's absent there,
            // or increments its count by 1 if it is already present.
            M.upsert(w, [](const uint64_t v)
            {
                return v + 1;
            }, 1, token);

            if(kmer_idx + k == len)
                break;

            w.advance(seq[kmer_idx + k]);   // Move the window to the next k-mer.
            kmer_idx++;
        }
    }

    std::cout << "#Unique k-mers: " << M.size() << "\n";
}


int main()
{
    hash_set_sequential();

    hash_set_concurrent();

    hash_map_sequential();

    return 0;
}
