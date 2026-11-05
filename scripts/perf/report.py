#!/usr/bin/env python3
"""Reports per-phase hardware-counter metrics from the runs of `run.sh` in a directory, as TSV lines of phase, metric,
and value. Per-op values divide by the phase's operation count; ratios of events use the run that counted both, else
the ratio of their per-op values. Runs whose hardware group was multiplexed are skipped, with a warning."""

import glob
import os
import sys

FIXED = ["cycles", "instructions", "task_clock_ns", "page_faults", "context_switches", "cpu_migrations"]

# Derived metrics: name, numerator, denominator; `None` denominators mean per op, and `cycles` gives cycle shares.
RATIOS = [
    ("mlp_l1d", "l1d_pend_miss.pending", "l1d_pend_miss.pending_cycles"),
    ("fb_full_share", "l1d_pend_miss.fb_full", "cycles"),
    ("mlp_dram", "offcore_requests_outstanding.demand_data_rd", "offcore_requests_outstanding.cycles_with_demand_data_rd"),
    ("l3_miss_latency_cycles", "offcore_requests_outstanding.l3_miss_demand_data_rd", "offcore_requests.l3_miss_demand_data_rd"),
    ("l3_stall_share", "cycle_activity.stalls_l3_miss", "cycles"),
    ("walk_cycles_per_walk", "dtlb_load_misses.walk_pending", "dtlb_load_misses.walk_completed"),
    ("walk_active_share", "dtlb_load_misses.walk_active", "cycles"),
    ("walk_2m_share", "dtlb_load_misses.walk_completed_2m_4m", "dtlb_load_misses.walk_completed"),
    ("sb_stall_share", "resource_stalls.sb", "cycles"),
]


def load(path):
    vals = {}
    with open(path) as f:
        for line in f:
            key, _, val = line.rstrip("\n").partition("\t")
            vals[key] = val
    return vals


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <out-dir>")

    runs = {os.path.basename(p)[:-4]: load(p) for p in sorted(glob.glob(os.path.join(sys.argv[1], "*.tsv")))}
    if not runs:
        sys.exit("No runs found.")

    first = next(iter(runs.values()))
    phases = [key[:-4] for key in first if key.endswith(".ops") and not key.startswith("perf.")]
    threads = float(first["threads"])

    def count(run, phase, event):
        val = run.get(f"perf.{phase}.{event}")
        return None if val is None or float(run.get(f"perf.{phase}.running_pct", 0)) < 99.9 else float(val)

    def per_op(phase, event):
        for run in runs.values():
            val = count(run, phase, event)
            if val is not None:
                return val / float(run[f"{phase}.ops"])
        return None

    def ratio(phase, num, den):
        for run in runs.values():
            n, d = count(run, phase, num), count(run, phase, den)
            if n is not None and d is not None:
                return n / d if d > 0 else float("nan")

        n, d = per_op(phase, num), per_op(phase, den)
        if n is not None and d is not None and den != "cycles":
            return n / d if d > 0 else float("nan")
        return None

    print("phase\tmetric\tvalue")
    for phase in phases:
        ops = float(first[f"{phase}.ops"])
        out = []

        valid = [run for run in runs.values() if count(run, phase, "cycles") is not None]
        for group, run in runs.items():
            if run not in valid:
                print(f"Skipping multiplexed group {group} at phase {phase}.", file=sys.stderr)
        if not valid:
            continue

        for event in FIXED:     # Counted in every run; averaged.
            out.append((f"{event}.per_op", sum(count(run, phase, event) for run in valid) / len(valid) / ops))

        out.append(("ipc", sum(count(r, phase, "instructions") / count(r, phase, "cycles") for r in valid) / len(valid)))
        out.append(("cpu_utilization", sum(count(r, phase, "task_clock_ns") / (threads * float(r[f"{phase}.seconds"]) * 1e9)
                                           for r in valid) / len(valid)))

        prefix = f"perf.{phase}."
        for run in valid:
            for key, val in run.items():
                event = key[len(prefix):]
                if key.startswith(prefix) and event not in FIXED and event != "running_pct":
                    out.append((f"{event}.per_op", float(val) / ops))

        for name, num, den in RATIOS:
            val = ratio(phase, num, den)
            if val is not None:
                out.append((name, val))

        demand, llc = per_op(phase, "mem_load_retired.l3_miss"), per_op(phase, "longest_lat_cache.miss")
        stalls = per_op(phase, "cycle_activity.stalls_l3_miss")
        if demand is not None and llc is not None:
            out.append(("noise_llc_miss.per_op", llc - demand))     # Prefetch, RFO, code, and page-walk misses.
        if demand is not None and stalls is not None:
            out.append(("stall_cycles_per_l3_miss", stalls / demand if demand > 0 else float("nan")))

        for name, val in out:
            print(f"{phase}\t{name}\t{val:.6g}")


if __name__ == "__main__":
    main()
