#!/usr/bin/env python3
"""Places runs of the benchmarks on a roofline, from the outputs of `scripts/bench/collate.py`: for each build and query
phase of the selected run configurations, its intensity, the operations per cache line that it fetches from beyond the
last-level cache (`longest_lat_cache.miss`); its rate; and the two roofs over it. The memory roof is the intensity
times `memory-roof`'s rate of random cache lines at the run's thread count: of updates for builds, as their lines are
written back too, and of reads for queries. The compute roof is the rate at which the run's threads, a core each, could
execute its instructions per operation, at its measured clock and the cores' issue width. Writes a row per phase, with
the roof that bounds it, the share of that roof that it reaches, and the shares of its cycles stalled, from runs that
counted the `stalls` group, if given; and prints them. Exits non-zero if a phase lacks its counters, or runs above a
roof by more than the tolerance, as the roof would then not bound it."""

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
FIELDS = CONFIG + ["phase", "roof", "mops", "lines_per_op", "ops_per_line", "line_mops", "roof_mops", "roof_share",
                   "memory_mops", "instructions_per_op", "ghz", "ipc", "cpu_utilization", "compute_mops", "ridge",
                   "bound", "efficiency"] + list(STALLS)


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
    per thread count; and the cache lines that each phase fetched per line that it accessed, per thread count and phase,
    where counted."""
    roofs = {}
    for row in read_tsv(os.path.join(d, "summary.tsv")):
        if row.get("random_read.mops") and row.get("random_update.mops"):
            roofs[int(row["threads"])] = {p: float(row[f"{p}.mops"]) for p in ROOF_PHASES if row.get(f"{p}.mops")}

    lines = {}
    for row in read_tsv(os.path.join(d, "counters.tsv")):
        if row["metric"] == LINES:
            lines[(int(row["threads"]), row["phase"])] = float(row["value"])

    return roofs, lines


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
                        help="instructions that a core can issue per cycle (default: 4)")
    parser.add_argument("--tolerance", type=float, default=0.05,
                        help="share by which a phase may run above a roof, as noise (default: 0.05)")
    parser.add_argument("-o", "--out", default="roofline.tsv", help="output file (default: roofline.tsv)")
    args = parser.parse_args()

    summaries, metrics = load_runs(args.dirs)
    roofs, roof_lines = load_roofs(args.roofs)
    stall_metrics = load_runs([args.stalls])[1] if args.stalls else {}
    t = args.threads
    if t not in roofs:
        sys.exit(f"No roofs at {t} threads in {args.roofs}")

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
            roof_mops = roofs[t][roof]
            lines = m[LINES]
            ghz = m["cycles_per_ref_cycle"] * args.base_ghz
            memory = roof_mops / lines
            compute = args.width * t * ghz * 1e3 / m["instructions.per_op"]
            row = dict(zip(CONFIG, cfg), phase=phase, roof=roof, mops=mops, lines_per_op=lines, ops_per_line=1 / lines,
                       line_mops=mops * lines, roof_mops=roof_mops, roof_share=mops * lines / roof_mops,
                       memory_mops=memory, instructions_per_op=m["instructions.per_op"], ghz=ghz, ipc=m.get("ipc"),
                       cpu_utilization=m.get("cpu_utilization"), compute_mops=compute, ridge=compute / roof_mops,
                       bound="memory" if memory <= compute else "compute", efficiency=mops / min(memory, compute))

            if mops > memory * (1 + args.tolerance):
                errors.append(f"{name} fetches {mops * lines:.6g} M lines/s, above its memory roof of {roof_mops:.6g}")
            if mops > compute * (1 + args.tolerance):
                errors.append(f"{name} runs at {mops:.6g} Mops, above its compute roof of {compute:.6g}")

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
    print(f"Roofs at {t} threads, in M lines/s: random reads {r['random_read']:.6g}, random updates "
          f"{r['random_update']:.6g}" + (f"; stream {r['stream']:.6g}" if "stream" in r else "") +
          (f", chase {r['chase']:.6g}" if "chase" in r else ""))
    head = f"{'table':<11} {'phase':<22} {'Mops':>7} {'ops/line':>8} {'roof':>5} {'memory':>7} {'compute':>7} " \
           f"{'bound':<7} {'effic.':>6} {'IPC':>5} {'busy':>5}" + (f" {'mem.stall':>9}" if args.stalls else "")
    print(head)
    pct = lambda v: f"{100 * v:.0f}%" if v is not None else "-"
    for row in rows:
        line = f"{row['table']:<11} {row['phase']:<22} {row['mops']:>7.0f} {row['ops_per_line']:>8.3f} " \
               f"{pct(row['roof_share']):>5} {row['memory_mops']:>7.0f} {row['compute_mops']:>7.0f} " \
               f"{row['bound']:<7} {pct(row['efficiency']):>6} {row['ipc'] or 0:>5.2f} {pct(row['cpu_utilization']):>5}"
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
