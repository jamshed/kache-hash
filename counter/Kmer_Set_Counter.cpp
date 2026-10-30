
#include "Kmer_Set_Counter.hpp"
#include "Fastx_Reader.hpp"
#include "Streaming_Kmer_Hash_Table.hpp"
#include "utility.hpp"

#include <stdexcept>
#include <thread>
#include <utility>


namespace kache_hash
{

// The counter of 31-mers, with minimizers of length 19.
struct Kmer_Set_Counter::Impl
{
    static constexpr uint16_t k = 31;
    static constexpr uint16_t l = 19;
    static constexpr std::size_t max_thread_c = 128;    // As many users as the table takes.

    typedef Streaming_Kmer_Hash_Table<k, true, uint32_t, l> table_t;

    const std::size_t thread_c;
    table_t M;
    std::vector<table_t::Token> token;  // The workers'.
    Timings timings;

    Impl(const std::size_t n, const std::size_t thread_c, const double lf):
          thread_c(thread_c)
        , M(n, thread_c, lf)
    {
        for(std::size_t t = 0; t < thread_c; ++t)
            token.push_back(M.register_user());
    }
};


bool Kmer_Set_Counter::supports(const uint16_t k)
{
    return k == Impl::k;
}


Kmer_Set_Counter::Kmer_Set_Counter(const uint16_t k, const std::size_t n, const std::size_t thread_c, const double lf)
{
    if(!supports(k))
        throw std::invalid_argument("k-mers of length " + std::to_string(k) + " are not supported, only of " + std::to_string(Impl::k));
    if(thread_c < 1 || thread_c > Impl::max_thread_c)
        throw std::invalid_argument("the thread count must be in [1, " + std::to_string(Impl::max_thread_c) + "]");

    const auto t_0 = now();
    impl = std::make_unique<Impl>(n, thread_c, lf);
    impl->timings.alloc = duration(now() - t_0);
}


void Kmer_Set_Counter::prime(const std::string& path)
{
    const auto t_0 = now();
    auto& M = impl->M;
    const auto& token = impl->token;
    Fastx_Reader(path).read(impl->thread_c, [&](const std::size_t t, const char* const* const seqs, const std::size_t* const lens, const std::size_t c)
    {
        M.upsert(seqs, lens, c, [](const uint32_t v){ return v; }, uint32_t(0), token[t]);
    });

    impl->timings.prime += duration(now() - t_0);
}


std::size_t Kmer_Set_Counter::count(const std::string& path)
{
    const auto t_0 = now();
    auto& M = impl->M;
    std::vector<Padded<std::size_t>> found(impl->thread_c, std::size_t(0));   // The workers'.
    Fastx_Reader(path).read(impl->thread_c, [&](const std::size_t t, const char* const* const seqs, const std::size_t* const lens, const std::size_t c)
    {
        found[t].unwrap() += M.update(seqs, lens, c, [](const uint32_t v){ return v + 1; });
    });

    std::size_t found_c = 0;
    for(auto& f : found)
        found_c += f.unwrap();

    impl->timings.count += duration(now() - t_0);
    return found_c;
}


std::vector<std::size_t> Kmer_Set_Counter::histogram(const std::size_t max)
{
    const auto t_0 = now();
    const auto& M = impl->M;
    const auto w_c = impl->thread_c;
    std::vector<std::vector<std::size_t>> H(w_c, std::vector<std::size_t>(max + 1, 0)); // The workers'.
    std::vector<std::thread> W;
    for(std::size_t w = 0; w < w_c; ++w)
        W.emplace_back([&, w]()
        {
            auto& h = H[w];
            M.for_each([&](const Kmer<Impl::k>&, const uint32_t v){ if(v > 0 && v <= max) h[v]++; }, w_c, w);
        });

    for(auto& w : W)
        w.join();

    for(std::size_t w = 1; w < w_c; ++w)
        for(std::size_t c = 1; c <= max; ++c)
            H[0][c] += H[w][c];

    impl->timings.histogram += duration(now() - t_0);
    return std::move(H[0]);
}


uint64_t Kmer_Set_Counter::count_of(const char* const kmer) const
{
    const auto v = impl->M.find(kmer);
    return v ? *v : 0;
}


std::size_t Kmer_Set_Counter::size() const
{
    return impl->M.size();
}


const Kmer_Set_Counter::Timings& Kmer_Set_Counter::timings() const
{
    return impl->timings;
}


}
