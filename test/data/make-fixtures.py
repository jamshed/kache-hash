#!/usr/bin/env python3
"""Makes the tests' small inputs in this directory from the datasets of `scripts/fetch-data.sh small`, deterministically:
- `ecoli.fa.gz`: E. coli K-12 MG1655's 3.90-4.20 Mbp, with four of its seven rRNA operons, its most repetitive window;
- `chr21.fa.gz`: two pieces of GRCh38's chr21, at 8.00-8.20 Mbp, a quarter of it an assembly gap of N's, and at
  12.50-12.60 Mbp, alpha satellite with 16% of its 31-mers distinct;
- `reads.fq.gz`: 2,000 reads simulated from both, of 100-250 bases, either strand, with 1% substitutions, some N's
  and some lowercase runs;
- `edge-cases.fa`: sequences that are empty, shorter than k, exactly k and k + 1 long, all N, N-split, lowercase,
  IUPAC-coded, low-complexity and multi-line.
Usage: test/data/make-fixtures.py <data-dir>"""

import gzip
import os
import random
import sys


def read_fasta(path):
    """Returns the sequence of the single-record FASTA file at `path`, upper-cased."""
    with gzip.open(path, "rt") as f:
        return "".join(line.strip() for line in f if not line.startswith(">")).upper()


def write_gz(path, text):
    """Writes `text` gzipped to `path`, with no timestamp, so that the bytes depend on the content alone."""
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, compresslevel=9) as f:
        f.write(text.encode())


def fasta(records, width=60):
    """Formats the `(name, sequence)` pairs `records` as FASTA, with lines of `width` bases."""
    return "".join(f">{name}\n" + "".join(s[i:i + width] + "\n" for i in range(0, len(s), width)) for name, s in records)


def revcomp(s):
    return s.translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))[::-1]


def simulate_reads(sources, n, rng):
    """Returns `n` reads of 100-250 bases from the sequences `sources`, as FASTQ."""
    out = []
    for i in range(n):
        s = sources[rng.randrange(len(sources))]
        length = rng.randint(100, 250)
        start = rng.randrange(len(s) - length)
        r = list(s[start:start + length])
        for j in range(length):
            u = rng.random()
            if u < 0.01:
                r[j] = rng.choice("ACGT")
            elif u < 0.012:
                r[j] = "N"
        r = "".join(r)
        if rng.random() < 0.5:
            r = revcomp(r)
        if rng.random() < 0.05:     # A lowercase run.
            b = rng.randrange(length // 2)
            r = r[:b] + r[b:b + 20].lower() + r[b + 20:]
        out.append(f"@read{i}\n{r}\n+\n{'I' * len(r)}\n")
    return "".join(out)


def edge_cases(rng):
    """Returns the edge cases' FASTA."""
    def bases(n):
        return "".join(rng.choice("ACGT") for _ in range(n))
    lower = bases(100)
    lower = lower[:30] + lower[30:70].lower() + lower[70:]
    iupac = list(bases(120))
    for j in range(5, len(iupac), 17):
        iupac[j] = rng.choice("RYKMSWBDHV")
    records = [
        ("empty", ""),
        ("shorter-than-k", bases(30)),
        ("exactly-k", bases(31)),
        ("k-plus-1", bases(32)),
        ("all-n", "N" * 100),
        ("n-split", "NNNNN" + bases(40) + "NN" + bases(35) + "N" + bases(31) + "NNNNNNNNNN" + bases(30)),
        ("lowercase", lower),
        ("all-lowercase", bases(80).lower()),
        ("iupac", "".join(iupac)),
        ("homopolymer", "A" * 200),
        ("dinucleotide", "CA" * 150),
        ("multi-line", bases(500)),
    ]
    return fasta(records)


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    data, here = sys.argv[1], os.path.dirname(os.path.abspath(__file__))
    rng = random.Random(31)

    ecoli = read_fasta(os.path.join(data, "ecoli.fa.gz"))[3_900_000:4_200_000]
    chr21 = read_fasta(os.path.join(data, "chr21.fa.gz"))
    chr21_gap, chr21_sat = chr21[8_000_000:8_200_000], chr21[12_500_000:12_600_000]
    write_gz(os.path.join(here, "ecoli.fa.gz"), fasta([("ecoli:3900000-4200000", ecoli)]))
    write_gz(os.path.join(here, "chr21.fa.gz"),
             fasta([("chr21:8000000-8200000", chr21_gap), ("chr21:12500000-12600000", chr21_sat)]))
    write_gz(os.path.join(here, "reads.fq.gz"),
             simulate_reads([ecoli, chr21_gap.replace("N", ""), chr21_sat], 2000, rng))
    with open(os.path.join(here, "edge-cases.fa"), "w") as f:
        f.write(edge_cases(rng))


if __name__ == "__main__":
    main()
