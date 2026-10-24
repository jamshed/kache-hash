#!/usr/bin/env python3
"""Reports per-phase operation statistics from the output of an instrumented driver run (`INSTRUMENT=ON`), as TSV lines
of phase, metric, and value. Per-op values divide by the phase's operation count. Run-union lines are the distinct
lines of each minimizer run, summed over the runs; they estimate the compulsory misses with reuse only within runs."""

import sys

LEVELS = ["primary", "secondary", "overflow"]
OUTCOMES = ["found", "placed", "absent"]
LINES = ["meta_lines", "slot_read_lines", "slot_write_lines", "overflow_lines"]
WINDOW = 13     # Number of k-mers sharing a minimizer occurrence at most: k - l + 1; longer runs are repeats.


def load(path):
    vals = {}
    with open(path) as f:
        for line in f:
            key, _, val = line.rstrip("\n").partition("\t")
            vals[key] = val
    return vals


def runs_by_length(vals, prefix):
    """Returns the number of runs and their distinct lines, keyed by run length, of the phase with key prefix `prefix`."""
    count, lines = {}, {}
    for key, val in vals.items():
        for name, out in (("runs.len.", count), ("runs.lines.", lines)):
            if key.startswith(prefix + name):
                out[int(key[len(prefix + name):].rstrip("+"))] = int(val)
    return count, lines


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <driver output TSV>")

    vals = load(sys.argv[1])
    if vals.get("instrument") != "1":
        sys.exit("Not an instrumented run.")

    num = lambda key: int(vals[key])
    phases = [key[:-4] for key in vals if key.endswith(".ops") and not key.startswith("inst.")]

    print("phase\tmetric\tvalue")
    for phase in phases:
        prefix, op = f"inst.{phase}.", ("update" if phase == "build" else "query")
        ops = num(prefix + op + ".ops")
        out = []

        for level in LEVELS:
            out.append((f"level.{level}", sum(num(f"{prefix}{op}.{level}.{o}") for o in OUTCOMES) / ops))
        for outcome in OUTCOMES:
            out.append((f"outcome.{outcome}", sum(num(f"{prefix}{op}.{l}.{outcome}") for l in LEVELS) / ops))

        for kind in LINES:
            out.append((f"{kind}.per_op", num(prefix + op + "." + kind) / ops))
        out.append(("lines.per_op", sum(num(prefix + op + "." + kind) for kind in LINES) / ops))
        for kind in LINES:
            out.append((f"run_{kind}.per_op", num(prefix + "runs." + kind) / ops))
        out.append(("run_lines.per_op", sum(num(prefix + "runs." + kind) for kind in LINES) / ops))

        compares = num(prefix + op + ".key_compares")
        visits = num(prefix + op + ".overflow_visits")
        out.append(("key_compares.per_op", compares / ops))
        out.append(("checksum_false_matches.per_compare", num(prefix + op + ".checksum_false_matches") / compares if compares else 0))
        out.append(("blind_reads.per_op", num(prefix + op + ".blind_reads") / ops))
        out.append(("overflow_probes.per_visit", num(prefix + op + ".overflow_probes") / visits if visits else 0))
        out.append(("overflow_locks.per_op", num(prefix + op + ".overflow_locks") / ops))

        run_c, len_sum = num(prefix + "runs.count"), num(prefix + "runs.len_sum")
        out.append(("run_len.mean", len_sum / run_c))
        out.append(("run_len.traffic_weighted", num(prefix + "runs.len_sq_sum") / len_sum))     # Mean run length of an op.
        count, lines = runs_by_length(vals, prefix)
        for r in range(1, WINDOW + 1):
            if r in count:
                out.append((f"run_lines.len_{r}", lines[r] / count[r]))     # Distinct lines per run of length `r`.
        long = [r for r in count if r > WINDOW]
        if long:
            out.append((f"run_lines.len_over_{WINDOW}", sum(lines[r] for r in long) / sum(count[r] for r in long)))
            out.append((f"ops_in_runs_over_{WINDOW}.share", sum(r * count[r] for r in long) / len_sum))

        if phase == "build":
            tries = num(prefix + "bucket_lock.tries")
            out.append(("bucket_lock.tries.per_op", tries / ops))
            out.append(("bucket_lock.fail_share", num(prefix + "bucket_lock.fails") / tries if tries else 0))
            out.append(("reader_lock.per_op", num(prefix + "reader_lock.acquisitions") / ops))
            out.append(("reader_lock.waits", num(prefix + "reader_lock.waits")))
            out.append(("reader_lock.retries", num(prefix + "reader_lock.retries")))
            out.append(("writer_lock.acquisitions", num(prefix + "writer_lock.acquisitions")))
            out.append(("writer_lock.wait_seconds", num(prefix + "writer_lock.wait_ns") / 1e9))
            out.append(("resize.count", num(prefix + "resize.count")))
            out.append(("resize.seconds", num(prefix + "resize.ns") / 1e9))
            placed = sum(num(prefix + "resize.placed." + l) for l in LEVELS)
            for level in LEVELS:
                out.append((f"resize.placed.{level}", num(prefix + "resize.placed." + level) / placed if placed else 0))
        out.append(("lines.unremembered", num(prefix + "lines.unremembered")))

        for name, val in out:
            print(f"{phase}\t{name}\t{val:.6g}")

    occupancy = {int(key.split(".")[-1]): int(val) for key, val in vals.items() if key.startswith("table.occupancy.")}
    if occupancy:
        buckets, full = sum(occupancy.values()), max(occupancy)
        print(f"table\tbuckets.empty_share\t{occupancy.get(0, 0) / buckets:.6g}")
        print(f"table\tbuckets.full_share\t{occupancy.get(full, 0) / buckets:.6g}")
        print(f"table\tbucket_size.mean\t{sum(s * c for s, c in occupancy.items()) / buckets:.6g}")
        print(f"table\toverflow_load\t{num('table.overflow_size') / num('table.overflow_capacity'):.6g}")


if __name__ == "__main__":
    main()
