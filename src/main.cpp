

#include "benchmark/benchmark.h"

void BM_file_read(benchmark::State&);

void scan_minimizer_correctness(benchmark::State&);
void BM_scan_minimizer_baseline(benchmark::State&);
void BM_scan_minimizer(benchmark::State&);
void BM_scan_minimizer_pos(benchmark::State&);

void BM_file_scan(benchmark::State&);
void BM_kmer_hash(benchmark::State&);
void BM_set_insert_generic_3rd_party(benchmark::State&);
void BM_set_insert_generic(benchmark::State&);
void BM_set_insert_batched(benchmark::State&);
void BM_set_insert_streaming(benchmark::State&);

template<bool track_stats_ = false> void BM_set_insert_streaming_blocked(benchmark::State&);


BENCHMARK(BM_file_read)->Unit(benchmark::kMillisecond);
// BENCHMARK(BM_scan_minimizer_baseline)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_scan_minimizer)->Unit(benchmark::kMillisecond);
// BENCHMARK(BM_scan_minimizer_pos)->Unit(benchmark::kMillisecond);

// BENCHMARK(BM_kmer_hash)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
// BENCHMARK(BM_set_insert_generic)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
// BENCHMARK(BM_set_insert_batched)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
// BENCHMARK(BM_set_insert_streaming)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);

BENCHMARK(BM_kmer_hash)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_set_insert_generic_3rd_party)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_set_insert_generic)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_set_insert_batched)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_set_insert_streaming)->Unit(benchmark::kMillisecond);

BENCHMARK(BM_set_insert_streaming_bucketed<false>)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_set_insert_streaming_bucketed<true>)->Unit(benchmark::kMillisecond);


BENCHMARK_MAIN();
