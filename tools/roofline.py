#!/usr/bin/env python3
"""Places runs of the benchmarks on an instruction roofline, from the outputs of `scripts/bench/collate.py`: for each
build and query phase of the selected run configurations, its performance, the instructions that its threads retire
per second; its intensity, the instructions that it executes per byte that it fetches from beyond the last-level cache
(`longest_lat_cache.miss`); and the roofs over it, the same for every table. The compute roof is the rate at which the
run's threads, a core each, could issue instructions, at the cores' issue width and the clock that `memory-roof`
measured at the run's thread count. The memory roofs are `memory-roof`'s rates at that thread count, in bytes: of
sequential reads, the hardware's limit; and of random cache lines, the limit of random access, updated for builds, as
their lines are written back too, and read for queries. The lower of the compute roof and the random-access roof bounds
a phase. Writes a row per phase, with the roof that bounds it, the share of that roof that it reaches, and the
shares of its cycles stalled, from runs that counted the `stalls` group, if given; and prints them. Exits non-zero if a
phase lacks its counters, or runs above a roof by more than the tolerance, as the roof would then not bound it."""

import argparse
import csv
import os
import re
import sys

CONFIG = ["workload", "table", "kind", "variant", "sizing", "threads"]
PHASE_MOPS = re.compile(r"(?P<phase>build|query\..+)\.mops")   # A phase's rate, as `summary.tsv` names it.
LINES = "longest_lat_cache.miss.per_op"     # Cache lines fetched from beyond the last-level cache, per operation.
ROOF_PHASES = ["random_read", "random_update", "stream", "chase"]  # The phases of `memory-roof`, in its output's order.
STALLS = {"stall_share": "cycle_activity.stalls_total.per_op",   # Shares of the cycles, by the stalls of the group.
          "mem_stall_share": "cycle_activity.stalls_mem_any.per_op",
          "l3_stall_share": "cycle_activity.stalls_l3_miss.per_op"}
FIELDS = CONFIG + ["phase", "roof", "mops", "instructions_per_op", "lines_per_op", "ops_per_line", "gips",
                   "instructions_per_byte", "gb_per_s", "roof_gb_per_s", "roof_share", "sequential_gb_per_s",
                   "sequential_share", "memory_gips", "compute_gips", "ridge", "bound", "efficiency", "ghz", "ipc",
                   "cpu_utilization"] + list(STALLS)


def read_tsv(path):
    """Returns the rows of the tab-separated file at `path`, as dictionaries keyed by its header."""
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def config(row):
    """Returns the run configuration of the collated row `row`, without its campaign."""
    return tuple(row[c] for c in CONFIG[:-1]) + (int(row["threads"]),)


def load_runs(dirs):
    """Returns the summary rows of the collated directories `dirs`, keyed by run configuration, and their counter
    metrics, keyed by run configuration and phase. Exits on a configuration in several of them."""
    summaries, metrics = {}, {}
    for d in dirs:
        for row in read_tsv(os.path.join(d, "summary.tsv")):
            cfg = config(row)
            if cfg in summaries:
                sys.exit(f"The configuration {'.'.join(map(str, cfg))} is in several directories")
            summaries[cfg] = row

        for row in read_tsv(os.path.join(d, "counters.tsv")):
            metrics.setdefault((config(row), row["phase"]), {})[row["metric"]] = float(row["value"])

    return summaries, metrics


def load_roofs(d):
    """Returns the rates of `memory-roof`'s phases in the collated directory `d`, in millions of cache lines per second,
    per thread count; the cache lines that each phase fetched per line that it accessed, per thread count and phase,
    where counted; and the highest clock of its phases, in cycles per reference cycle, per thread count, where
    counted."""
    roofs = {}
    for row in read_tsv(os.path.join(d, "summary.tsv")):
        if row.get("random_read.mops") and row.get("random_update.mops"):
            roofs[int(row["threads"])] = {p: float(row[f"{p}.mops"]) for p in ROOF_PHASES if row.get(f"{p}.mops")}

    lines, clocks = {}, {}
    for row in read_tsv(os.path.join(d, "counters.tsv")):
        t = int(row["threads"])
        if row["metric"] == LINES:
            lines[(t, row["phase"])] = float(row["value"])
        elif row["metric"] == "cycles_per_ref_cycle" and row["phase"] in ROOF_PHASES:
            clocks[t] = max(clocks.get(t, 0), float(row["value"]))

    return roofs, lines, clocks


def phases(summary):
    """Returns the build and query phases of the summary row `summary` that ran at its thread count, with their rates;
    not the queries at further thread counts, `query.<input>.T<count>`, whose roofs differ."""
    out = []
    for col, val in summary.items():
        m = PHASE_MOPS.fullmatch(col or "")
        if m and val and not re.search(r"\.T\d+$", m["phase"]):
            out.append((m["phase"], float(val)))

    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split(":")[0] + ".")
    parser.add_argument("dirs", nargs="+", help="collated directories of the benchmarks' runs")
    parser.add_argument("--roofs", required=True, help="collated directory of `memory-roof`'s runs")
    parser.add_argument("--stalls", help="collated directory of runs that counted the `stalls` group")
    parser.add_argument("--workload", default="reference", help="workload of the runs (default: reference)")
    parser.add_argument("--kind", default="set", help="kind of the tables (default: set)")
    parser.add_argument("--variant", default="mt", help="variant of the tables (default: mt)")
    parser.add_argument("--sizing", default="presized",
                        help="sizing of the tables, besides `na` for those that can not be sized (default: presized)")
    parser.add_argument("--threads", type=int, default=56, help="thread count of the runs (default: 56)")
    parser.add_argument("--base-ghz", type=float, default=2.2,
                        help="reference clock, at which `ref_cycles` counts, in GHz (default: 2.2)")
    parser.add_argument("--width", type=int, default=4,
                        help="instructions that a core can issue per cycle, for the compute roof (default: 4)")
    parser.add_argument("--line-bytes", type=int, default=64, help="bytes in a cache line (default: 64)")
    parser.add_argument("--tolerance", type=float, default=0.05,
                        help="share by which a phase may run above a roof, as noise (default: 0.05)")
    parser.add_argument("-o", "--out", default="roofline.tsv", help="output file (default: roofline.tsv)")
    args = parser.parse_args()

    summaries, metrics = load_runs(args.dirs)
    roofs, roof_lines, clocks = load_roofs(args.roofs)
    stall_metrics = load_runs([args.stalls])[1] if args.stalls else {}
    t = args.threads
    if t not in roofs or "stream" not in roofs[t]:
        sys.exit(f"No roofs at {t} threads in {args.roofs}")
    if t not in clocks:
        sys.exit(f"No clock of `memory-roof` at {t} threads in {args.roofs}, for the compute roof")

    gb = lambda line_mops: line_mops * args.line_bytes / 1e3    # M cache lines/s, in GB/s.
    roof_ghz = clocks[t] * args.base_ghz
    compute = args.width * t * roof_ghz     # The compute roof, in G instructions/s.
    sequential = gb(roofs[t]["stream"])

    errors, warnings, rows = [], [], []
    for p in ("random_read", "random_update"):
        n = roof_lines.get((t, p))
        if n is None or abs(n - 1) > 0.1:
            warnings.append(f"`memory-roof`'s {p} fetched {n if n is not None else 'uncounted'} lines per line"
                            f" accessed, not about 1: its rate and the intensities count different lines")

    for cfg, summary in sorted(summaries.items()):
        workload, table, kind, variant, sizing, threads = cfg
        if (workload, kind, variant, threads) != (args.workload, args.kind, args.variant, t) or \
                sizing not in (args.sizing, "na"):
            continue

        for phase, mops in sorted(phases(summary)):
            name = f"{table} {phase}"
            m = metrics.get((cfg, phase), {})
            missing = [x for x in (LINES, "instructions.per_op", "cycles_per_ref_cycle") if x not in m]
            if missing:
                errors.append(f"{name} lacks the counters {', '.join(missing)}")
                continue

            roof = "random_update" if phase == "build" else "random_read"
            roof_gb = gb(roofs[t][roof])
            lines, instructions = m[LINES], m["instructions.per_op"]
            gips, gb_per_s = mops * instructions / 1e3, gb(mops * lines)
            per_byte = instructions / (lines * args.line_bytes)
            memory = roof_gb * per_byte
            row = dict(zip(CONFIG, cfg), phase=phase, roof=roof, mops=mops, instructions_per_op=instructions,
                       lines_per_op=lines, ops_per_line=1 / lines, gips=gips, instructions_per_byte=per_byte,
                       gb_per_s=gb_per_s, roof_gb_per_s=roof_gb, roof_share=gb_per_s / roof_gb,
                       sequential_gb_per_s=sequential, sequential_share=gb_per_s / sequential, memory_gips=memory,
                       compute_gips=compute, ridge=compute / roof_gb,
                       bound="memory" if memory <= compute else "compute",
                       efficiency=gips / min(memory, compute), ghz=m["cycles_per_ref_cycle"] * args.base_ghz,
                       ipc=m.get("ipc"), cpu_utilization=m.get("cpu_utilization"))

            if gb_per_s > roof_gb * (1 + args.tolerance):
                errors.append(f"{name} fetches {gb_per_s:.6g} GB/s, above its memory roof of {roof_gb:.6g}")
            if gips > compute * (1 + args.tolerance):
                errors.append(f"{name} runs {gips:.6g} G instructions/s, above the compute roof of {compute:.6g}")

            if args.stalls:
                s = stall_metrics.get((cfg, phase), {})
                if "cycles.per_op" not in s:
                    warnings.append(f"{name} has no stall counts in {args.stalls}")
                for col, event in STALLS.items():
                    if event in s and "cycles.per_op" in s:
                        row[col] = s[event] / s["cycles.per_op"]

            rows.append(row)

    with open(args.out, "w") as f:
        print("\t".join(FIELDS), file=f)
        for row in rows:
            cells = [row.get(c) for c in FIELDS]
            print("\t".join("" if v is None else f"{v:.6g}" if isinstance(v, float) else str(v) for v in cells), file=f)

    r = roofs[t]
    print(f"Roofs at {t} threads: compute {compute:.6g} G instructions/s ({args.width} per cycle at {roof_ghz:.2f}"
          f" GHz); memory, in GB/s: sequential reads {sequential:.6g}, random reads {gb(r['random_read']):.6g}, random"
          f" updates {gb(r['random_update']):.6g}" + (f", chase {gb(r['chase']):.6g}" if "chase" in r else ""))
    head = f"{'table':<16} {'phase':<22} {'Mops':>7} {'GIPS':>7} {'instr/B':>7} {'GB/s':>6} {'roof':>5} " \
           f"{'memory':>7} {'bound':<7} {'effic.':>6} {'IPC':>5} {'busy':>5}" + \
           (f" {'mem.stall':>9}" if args.stalls else "")
    print(head)
    pct = lambda v: f"{100 * v:.0f}%" if v is not None else "-"
    for row in rows:
        line = f"{row['table']:<16} {row['phase']:<22} {row['mops']:>7.0f} {row['gips']:>7.1f} " \
               f"{row['instructions_per_byte']:>7.2f} {row['gb_per_s']:>6.1f} {pct(row['roof_share']):>5} " \
               f"{row['memory_gips']:>7.0f} {row['bound']:<7} {pct(row['efficiency']):>6} {row['ipc'] or 0:>5.2f} " \
               f"{pct(row['cpu_utilization']):>5}"
        if args.stalls:
            line += f" {pct(row.get('mem_stall_share')):>9}"
        print(line)

    print(f"Wrote {len(rows)} phases to {args.out}.")
    for w in warnings:
        print(f"Warning: {w}.", file=sys.stderr)
    for e in errors:
        print(f"Error: {e}.", file=sys.stderr)

    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
