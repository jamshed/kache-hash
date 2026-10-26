#!/usr/bin/env python3
"""Collates the runs of `run.sh` in an output directory into `summary.tsv`, a row per workload and run configuration
with the medians over its repetitions, and `build-cost.tsv`, SSHash's index pipeline (GGCAT, then SSHash) against
the dynamic tables' load, init, and build at the same thread count. Exits non-zero if the runs of a workload disagree
on the table size or on any query's hits, or if a build's operation count is not its input's k-mer count."""

import glob
import os
import re
import statistics
import sys

RUN = re.compile(r"(?P<table>[^.]+)\.(?P<kind>set|map)\.(?P<variant>mt|seq)\.(?P<sizing>[^.]+)\.T(?P<threads>\d+)"
                 r"\.r(?P<rep>\d+)\.tsv")
CONFIG = ["workload", "table", "kind", "variant", "sizing", "threads"]

# Measurements reported as medians: per-run keys, and per-query key suffixes.
MEDIANS = ["load.seconds", "init.seconds", "build.seconds", "build.mops", "memory.build_max_rss_bytes"]
QUERY_MEDIANS = ["seconds", "mops"]


def load(path):
    vals = {}
    with open(path) as f:
        for line in f:
            key, _, val = line.rstrip("\n").partition("\t")
            vals[key] = val
    return vals


def median(runs, key):
    vals = [float(run[key]) for run in runs if key in run]
    if not vals:
        return ""
    m = statistics.median(vals)
    return str(int(m)) if m.is_integer() else f"{m:.6g}"


def name(cfg):
    """Returns the name of the run configuration `cfg` within its workload, as `run.sh` names its outputs."""
    return f"{cfg[1]}.{cfg[2]}.{cfg[3]}.{cfg[4]}.T{cfg[5]}"


def where(cfgs):
    """Returns a description of the run configurations `cfgs`: their names, if few."""
    return ", ".join(sorted(cfgs)) if len(cfgs) <= 3 else f"{len(cfgs)} configurations"


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <out-dir>")

    out_dir = sys.argv[1]
    configs = {}    # Runs per configuration.
    costs = {}  # SSHash index costs per workload.
    for wl_dir in sorted(glob.glob(os.path.join(out_dir, "*", ""))):
        wl = os.path.basename(os.path.dirname(wl_dir))
        for path in sorted(glob.glob(os.path.join(wl_dir, "*.tsv"))):
            m = RUN.fullmatch(os.path.basename(path))
            if m:
                cfg = (wl, m["table"], m["kind"], m["variant"], m["sizing"], int(m["threads"]))
                configs.setdefault(cfg, []).append(load(path))
            elif os.path.basename(path) == "index.cost.tsv":
                cost = load(path)
                if "sshash.index_bytes" in cost:    # Its last line; a failed build leaves it out.
                    costs[wl] = cost

    if not configs:
        sys.exit("No runs found.")

    errors = []
    for wl in sorted({cfg[0] for cfg in configs}):
        runs = [(cfg, run) for cfg, R in configs.items() if cfg[0] == wl for run in R]
        keys = {"table.size"} | {key for _, run in runs for key in run if re.fullmatch(r"query\d+\.hits", key)}
        for key in sorted(keys):
            vals = {}
            for cfg, run in runs:
                vals.setdefault(run.get(key, "missing"), set()).add(name(cfg))
            if len(vals) > 1:
                errors.append(f"{wl}: {key} differs: " + "; ".join(f"{v} in {where(c)}" for v, c in vals.items()))

        for cfg, run in runs:
            if "build.ops" in run and run["build.ops"] != run["index.kmers"]:
                errors.append(f"{wl}: {name(cfg)} built {run['build.ops']} of {run['index.kmers']} k-mers")

    query_c = max(sum(1 for key in run if re.fullmatch(r"query\d+\.ops", key))
                  for R in configs.values() for run in R)
    query_cols = [f"query{q}.{s}" for q in range(query_c) for s in QUERY_MEDIANS]
    with open(os.path.join(out_dir, "summary.tsv"), "w") as f:
        print("\t".join(CONFIG + ["reps", "table.size"] + MEDIANS + query_cols), file=f)
        for cfg, R in sorted(configs.items()):
            cells = [str(c) for c in cfg] + [str(len(R)), R[0].get("table.size", "")]
            print("\t".join(cells + [median(R, key) for key in MEDIANS + query_cols]), file=f)

    with open(os.path.join(out_dir, "build-cost.tsv"), "w") as f:
        print("workload\tstructure\tsizing\tthreads\tseconds\tmax_rss_bytes\tstages", file=f)
        for wl, cost in sorted(costs.items()):
            gg, ss = float(cost["ggcat.seconds"]), float(cost["sshash.seconds"])
            rss = max(int(cost["ggcat.max_rss_bytes"]), int(cost["sshash.max_rss_bytes"]))
            print(f"{wl}\tsshash\tna\t{cost['threads']}\t{gg + ss:.6g}\t{rss}\tggcat {gg:.6g} + sshash {ss:.6g}", file=f)
            for cfg, R in sorted(configs.items()):
                if cfg[0] == wl and cfg[1] != "sshash" and cfg[2] == "set" and cfg[5] == int(cost["threads"]):
                    ld, it, bd = (float(median(R, f"{p}.seconds")) for p in ("load", "init", "build"))
                    rss = median(R, "memory.build_max_rss_bytes")
                    print(f"{wl}\t{cfg[1]}.{cfg[3]}\t{cfg[4]}\t{cfg[5]}\t{ld + it + bd:.6g}\t{rss}"
                          f"\tload {ld:.6g} + init {it:.6g} + build {bd:.6g}", file=f)

    print(f"Collated {sum(len(R) for R in configs.values())} runs of {len(configs)} configurations into {out_dir}.")
    for error in errors:
        print(error, file=sys.stderr)
    if errors:
        sys.exit(f"{len(errors)} disagreements.")


if __name__ == "__main__":
    main()
