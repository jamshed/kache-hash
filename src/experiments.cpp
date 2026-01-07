
#include "Streaming_Kmer_Hash_Table.hpp"
#include "Parser.hpp"

#include <cstdint>
#include <cstddef>
#include <string>

template <uint16_t k>
void throughput_kache_hash_sequential(const std::string& path, const std::size_t final_sz)
{
    Streaming_Kmer_Hash_Table<k, false, uint64_t> M(final_sz, 1);

    const auto t_s = now();

    uint64_t n = 0;
#ifndef QUIET
    constexpr uint64_t log_step = 32768;
    uint64_t checkpoint = log_step;
    std::cerr << "\n";
#endif

    Parser input(path);
    const auto token = M.register_user();
    while(input.read_next_seq())
    {
        auto const seq = input.seq();
        auto const len = input.seq_len();
        if(len < k)
            continue;

        Kmer_Window<k> w;
        w.init(seq);
        std::size_t kmer_idx = 0;
        while(true)
        {
            M.upsert(w, [](const uint64_t v)
            {
                return v + 1;
            }, 1, token);

            if(kmer_idx + k == len)
                break;

            w.advance(seq[kmer_idx + k]);
            kmer_idx++;
        }

        n += len - k + 1;
#ifndef QUIET
        if(n >= checkpoint)
        {
            std::cerr << "\rInserted " << n << " k-mers.";
            checkpoint = n + log_step;
        }
#endif
    }

    input.close();
    const auto t_e = now();
    std::cerr << "\n Time taken for kache-hash-sequential: " << duration(t_e - t_s) << "s.\n";

    return M.size();
}


int main(int argc, char* argv[])
{
    (void)argc, (void)argv;

    return 0;
}
