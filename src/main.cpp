

#include "benchmark/benchmark.h"

void BM_file_read(benchmark::State&);

void scan_minimizer_correctness(benchmark::State&);
void BM_scan_minimizer_baseline(benchmark::State&);
void BM_scan_minimizer(benchmark::State&);
void BM_scan_minimizer_pos(benchmark::State&);

void BM_file_scan(benchmark::State&);
void BM_kmer_hash(benchmark::State&);
void BM_set_insert_generic(benchmark::State&);
void BM_set_insert_batched(benchmark::State&);
void BM_set_insert_streaming(benchmark::State&);

template<bool track_stats_ = false> void BM_set_insert_streaming_blocked(benchmark::State&);


// BENCHMARK(BM_file_read)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);

// BENCHMARK(scan_minimizer_correctness);
// BENCHMARK(BM_scan_minimizer_baseline)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
// BENCHMARK(BM_scan_minimizer)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
// BENCHMARK(BM_scan_minimizer_pos)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);

// BENCHMARK(BM_kmer_hash)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
// BENCHMARK(BM_set_insert_generic)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
// BENCHMARK(BM_set_insert_batched)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
// BENCHMARK(BM_set_insert_streaming)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);

BENCHMARK(BM_file_scan)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_scan_minimizer)->Unit(benchmark::kMillisecond);

BENCHMARK(BM_set_insert_streaming_blocked<true>)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_set_insert_streaming_blocked<false>)->MinWarmUpTime(1)->MinTime(2)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
