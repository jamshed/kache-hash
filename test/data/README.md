# Test inputs

The small inputs of the `quick` tests (see `../CMakeLists.txt`), 248 KB in all. `make-fixtures.py` makes them,
deterministically, from the datasets of `scripts/fetch-data.sh small`: `make-fixtures.py <data-dir>` rewrites them
byte for byte.

| File | Content | Source |
|---|---|---|
| `ecoli.fa.gz` | 3.90–4.20 Mbp of E. coli K-12 MG1655, its most repetitive window, with four of its seven rRNA operons | NCBI RefSeq GCF_000005845.2 |
| `chr21.fa.gz` | chr21 at 8.00–8.20 Mbp, a quarter of it an assembly gap of N's, and at 12.50–12.60 Mbp, alpha satellite with 16% of its 31-mers distinct | GRCh38, UCSC hg38 |
| `reads.fq.gz` | 2,000 reads of 100–250 bases simulated from both, from either strand, with 1% substitutions, some N's and some lowercase runs | Simulated |
| `edge-cases.fa` | Sequences that are empty, shorter than k, exactly k and k + 1 long, all N, N-split, lowercase, IUPAC-coded, low-complexity and multi-line | Synthetic |
