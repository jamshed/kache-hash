

#include "benchmark/benchmark.h"

void BM_file_read(benchmark::State&);
void BM_kmer_hash(benchmark::State&);
void BM_set_insert_generic(benchmark::State&);
void BM_set_insert_batched(benchmark::State&);

void BM_scan_minimizer_baseline(benchmark::State& state);
void BM_scan_minimizer(benchmark::State& state);


BENCHMARK(BM_file_read)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_kmer_hash)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_set_insert_generic)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_set_insert_batched)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);

BENCHMARK(BM_scan_minimizer_baseline)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_scan_minimizer)->MinWarmUpTime(1)->MinTime(5)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
