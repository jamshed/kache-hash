

#include "benchmark/benchmark.h"

void BM_set_insert_generic(benchmark::State&);

BENCHMARK(BM_set_insert_generic)->MinWarmUpTime(1)->MinTime(5);

BENCHMARK_MAIN();
