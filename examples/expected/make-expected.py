#!/usr/bin/env python3
"""Makes the examples' expected outputs in this directory from the inputs in `test/data/`, by a plain count of their
31-mers, independent of kache-hash, in the examples' formats:
- `count_kmers.<input>.txt`: the distinct canonical 31-mers of `ecoli.fa.gz`, `reads.fq.gz` and `edge-cases.fa`, their
  occurrences, and the 31-mers by their counts, from 1 to 9 and 10 or more;
- `filter_reads.names`: the names of the reads of `reads.fq.gz` with at least half of their 31-mer positions in
  `ecoli.fa.gz`, in either strand;
- `strand_membership.tsv`: each read's 31-mer positions, and its 31-mers in `ecoli.fa.gz`'s forward and reverse strands;
- `count_in_set.txt`: the distinct canonical 31-mers of `ecoli.fa.gz`, the occurrences of the reads' 31-mers among them,
  and the former by their counts in the reads, from 1 to 10.
Non-ACGT characters split the sequences, and lowercase bases count as bases.
Usage: examples/expected/make-expected.py"""

import gzip
import os
from collections import Counter

K = 31
COMPLEMENT = str.maketrans("ACGT", "TGCA")
ACGT = str.maketrans("", "", "ACGT")


def records(path):
    """Returns the `(name, sequence)` records of the FASTA or FASTQ file at `path`, gzipped or not, as the examples read
    them: a FASTQ record has four lines, and a FASTA record's sequence joins its lines."""
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as f:
        lines = [line.rstrip("\r") for line in f.read().split("\n")]
    if lines[0].startswith("@"):
        return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4) if lines[i].startswith("@")]

    recs = []
    for line in lines:
        if line.startswith(">"):
            recs.append((line[1:], []))
        elif recs:
            recs[-1][1].append(line)
    return [(name, "".join(seq)) for name, seq in recs]


def kmers(seq):
    """Yields the 31-mers of the sequence `seq`, upper-cased, but those with a non-ACGT character."""
    seq = seq.upper()
    for i in range(len(seq) - K + 1):
        w = seq[i:i + K]
        if not w.translate(ACGT):
            yield w


def revcomp(w):
    return w.translate(COMPLEMENT)[::-1]


def canonical(w):
    return min(w, revcomp(w))


def positions(seq):
    """Returns the number of 31-mer positions of the sequence `seq`, with non-ACGT characters or not."""
    return max(len(seq) - K + 1, 0)


def count_kmers(path):
    counts = Counter(canonical(w) for _, s in records(path) for w in kmers(s))
    hist = Counter(min(c, 10) for c in counts.values())
    return ([f"distinct {K}-mers\t{len(counts)}", f"occurrences\t{sum(counts.values())}"] +
            [f"count {c}{'+' if c == 10 else ''}\t{hist[c]}" for c in range(1, 11)])


def filter_reads(ref, reads, threshold):
    present = {canonical(w) for _, s in records(ref) for w in kmers(s)}
    return [name for name, s in records(reads)
            if positions(s) > 0 and sum(canonical(w) in present for w in kmers(s)) >= threshold * positions(s)]


def strand_membership(ref, reads):
    # The reference's 31-mers as they occur in it: those of its forward strand. A read's 31-mer is in the reverse
    # strand iff its reverse complement is in the forward one.
    forward = {w for _, s in records(ref) for w in kmers(s)}
    lines = ["read\tpositions\tforward\treverse"]
    for i, (_, s) in enumerate(records(reads)):
        ws = list(kmers(s))
        lines.append(f"{i}\t{positions(s)}\t{sum(w in forward for w in ws)}\t{sum(revcomp(w) in forward for w in ws)}")
    return lines


def count_in_set(ref, reads):
    counts = {canonical(w): 0 for _, s in records(ref) for w in kmers(s)}
    counted = 0
    for _, s in records(reads):
        for w in kmers(s):
            c = canonical(w)
            if c in counts:
                counts[c] += 1
                counted += 1
    hist = Counter(counts.values())
    return ([f"set {K}-mers\t{len(counts)}", f"occurrences counted\t{counted}"] +
            [f"count {c}\t{hist[c]}" for c in range(1, 11)])


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    data = os.path.join(here, "..", "..", "test", "data")
    ecoli, reads, edges = (os.path.join(data, f) for f in ("ecoli.fa.gz", "reads.fq.gz", "edge-cases.fa"))
    outputs = {
        "count_kmers.ecoli.txt": count_kmers(ecoli),
        "count_kmers.reads.txt": count_kmers(reads),
        "count_kmers.edge-cases.txt": count_kmers(edges),
        "filter_reads.names": filter_reads(ecoli, reads, 0.5),
        "strand_membership.tsv": strand_membership(ecoli, reads),
        "count_in_set.txt": count_in_set(ecoli, reads),
    }
    for name, lines in outputs.items():
        with open(os.path.join(here, name), "w") as f:
            f.write("".join(line + "\n" for line in lines))


if __name__ == "__main__":
    main()
