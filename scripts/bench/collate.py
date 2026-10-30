#!/usr/bin/env python3
"""Collates the runs of `run.sh` in output directories, a campaign each, named after its directory, into the first
one's `summary.tsv`, a row per campaign, workload, and run configuration with the versions of the code that its runs
ran and the medians over its repetitions; `build-cost.tsv`, SSHash's index pipeline (GGCAT, then SSHash) against the
dynamic tables' load, init, and build at the same thread count; and `counters.tsv`, a row per campaign, workload, run
configuration, phase, and hardware-counter metric, with its median over the runs that counted it. Queries are named
after their inputs, as `query.<input>`. Exits non-zero if the runs of a workload, across the campaigns, disagree on
the table size or on any query's hits, if a query's hits differ across its thread counts, or if a build's operation
count is not its input's k-mer count. Warns of campaigns whose runs ran several versions of the code, of runs left with
much of their memory in small pages after the build, as on a node short of free huge pages, and of multiplexed
counters, whose counts are left out."""

import glob
import os
import re
import statistics
import sys

multi = False   # Whether there are several campaigns.

RUN = re.compile(r"(?P<table>[^.]+)\.(?P<kind>set|map)\.(?P<variant>mt|seq)\.(?P<sizing>[^.]+)\.T(?P<threads>\d+)"
                 r"\.r(?P<rep>\d+)\.tsv")
CONFIG = ["campaign", "workload", "table", "kind", "variant", "sizing", "threads"]
# A query phase's measurement, as the drivers name it, and an input's file extensions.
QUERY = re.compile(r"query(?P<q>\d+)(?P<threads>\.T\d+)?\.(?P<stat>hits|seconds|mops)")
INPUT_EXT = re.compile(r"(\.(fa|fna|fasta|fastq|fq))?(\.gz)?$")
# A phase that counts hardware events, as the drivers name it: the build, or a query at the run's thread count or a
# further one.
PHASE = re.compile(r"build|query(?P<q>\d+)(?P<threads>\.T(?P<t>\d+))?")

# Measurements reported as medians: per-run keys, and per-query key suffixes.
MEDIANS = ["load.seconds", "init.seconds", "build.seconds", "build.mops", "memory.build_max_rss_bytes",
           "memory.build_huge_share"]
QUERY_MEDIANS = ["seconds", "mops"]

# Runs warned of for small pages: those with at least 1 GiB of anonymous memory after the build, and less than 80% of
# it in transparent huge pages.
HUGE_MIN_KB = 1 << 20
HUGE_MIN_SHARE = 0.8


def load(path):
    vals = {}
    with open(path) as f:
        for line in f:
            key, _, val = line.rstrip("\n").partition("\t")
            vals[key] = val
    return vals


def name_queries(run):
    """Names the query phases of the run `run` after their inputs, keying the measurements of `query<i>` by
    `query.<input>`, the input's file name without its extensions. Returns the inputs whose hits differ between their
    phases at the run's thread count and at further ones (`query<i>.T<t>`)."""
    hits = {}   # Hits of each input's phases.
    for key in list(run):
        m = QUERY.fullmatch(key)
        if m and f"query{m['q']}" in run:
            name = INPUT_EXT.sub("", os.path.basename(run[f"query{m['q']}"]))
            if m["stat"] == "hits":
                hits.setdefault(name, set()).add(run[key])
            if not m["threads"]:
                run[f"query.{name}.{m['stat']}"] = run[key]

    return sorted(name for name, H in hits.items() if len(H) > 1)


def counters(run):
    """Returns the hardware-counter metrics of the run `run` per phase, named as in the summary: each event per
    operation, instructions per cycle, cycles per reference cycle (the effective frequency over the base one), and the
    share of the time that its threads spent on a CPU. Returns also the phases whose counters were multiplexed, which
    are left out."""
    phases = [key[len("perf."):-len(".running_pct")] for key in run
              if key.startswith("perf.") and key.endswith(".running_pct")]
    metrics, multiplexed = {}, []
    for phase in phases:
        m = PHASE.fullmatch(phase)
        if not m:
            continue

        name, threads = phase, float(run["threads"])
        if m["q"]:
            name = "query." + INPUT_EXT.sub("", os.path.basename(run[f"query{m['q']}"])) + (m["threads"] or "")
            threads = float(m["t"] or threads)
        if float(run[f"perf.{phase}.running_pct"]) < 99.9:
            multiplexed.append(name)
            continue

        events = {}     # The phase's counts, less those of the phases it prefixes, e.g. `query0.T1` of `query0`.
        for key, val in run.items():
            owner = max((p for p in phases if key.startswith(f"perf.{p}.")), key=len, default=None)
            if owner == phase and key != f"perf.{phase}.running_pct":
                events[key[len(f"perf.{phase}."):]] = float(val)

        ops = float(run[f"{phase}.ops"])
        vals = {f"{event}.per_op": count / ops for event, count in events.items()}
        if events.get("cycles"):
            vals["ipc"] = events["instructions"] / events["cycles"]
            if events.get("ref_cycles"):
                vals["cycles_per_ref_cycle"] = events["cycles"] / events["ref_cycles"]
        if "task_clock_ns" in events:
            vals["cpu_utilization"] = events["task_clock_ns"] / (threads * float(run[f"{phase}.seconds"]) * 1e9)
        metrics[name] = vals

    return metrics, multiplexed


def median(runs, key):
    vals = [float(run[key]) for run in runs if key in run]
    if not vals:
        return ""
    m = statistics.median(vals)
    return str(int(m)) if m.is_integer() else f"{m:.6g}"


def name(cfg):
    """Returns the name of the run configuration `cfg` within its workload, as `run.sh` names its outputs, with its
    campaign if there are several."""
    return f"{cfg[1]}.{cfg[2]}.{cfg[3]}.{cfg[4]}.T{cfg[5]}" + (f" ({cfg[6]})" if multi else "")


def known(cfg):
    """Returns whether the run configuration `cfg` has a known defect, that its disagreements are reported for but do
    not fail: iceberg's resizes duplicate and lose keys, even with one user thread, and its presized concurrent builds
    rarely hold a key that queries then miss."""
    return cfg[1] == "iceberg"


def where(cfgs):
    """Returns a description of the run configurations `cfgs`: their names, if few."""
    return ", ".join(sorted(cfgs)) if len(cfgs) <= 3 else f"{len(cfgs)} configurations"


def main():
    global multi
    if len(sys.argv) < 2:
        sys.exit(f"Usage: {sys.argv[0]} <out-dir> [<out-dir>...]")

    out_dirs = sys.argv[1:]
    multi = len(out_dirs) > 1
    configs = {}    # Runs per configuration.
    costs = {}  # SSHash index costs per workload and campaign.
    versions = {}   # Run counts per campaign and version.
    counted = {}    # Counter metrics of the runs per configuration, for the runs that counted any.
    errors, defects, warnings = [], [], []
    for out_dir in out_dirs:
        campaign = os.path.basename(os.path.normpath(out_dir))
        for wl_dir in sorted(glob.glob(os.path.join(out_dir, "*", ""))):
            wl = os.path.basename(os.path.dirname(wl_dir))
            for path in sorted(glob.glob(os.path.join(wl_dir, "*.tsv"))):
                m = RUN.fullmatch(os.path.basename(path))
                if m:
                    cfg = (wl, m["table"], m["kind"], m["variant"], m["sizing"], int(m["threads"]), campaign)
                    run = load(path)
                    for query in name_queries(run):
                        errors.append(f"{wl}: {query}'s hits differ across thread counts in {os.path.basename(path)}")

                    counts = versions.setdefault(campaign, {})
                    counts[run.get("version", "unknown")] = counts.get(run.get("version", "unknown"), 0) + 1

                    anon = int(run.get("memory.build_anon_kb", 0))
                    if anon > 0:    # Runs on Linux report the huge pages of their memory.
                        share = int(run.get("memory.build_anon_huge_pages_kb", 0)) / anon
                        run["memory.build_huge_share"] = f"{share:.3f}"
                        if anon >= HUGE_MIN_KB and share < HUGE_MIN_SHARE:
                            warnings.append(f"{wl}: huge pages back {share:.0%} of the {anon / (1 << 20):.1f} GiB of"
                                            f" anonymous memory after the build in {os.path.basename(path)}"
                                            + (f" ({campaign})" if multi else ""))

                    metrics, multiplexed = counters(run)
                    for phase in multiplexed:
                        warnings.append(f"{wl}: counters multiplexed in {phase} of {os.path.basename(path)}, so left"
                                        " out" + (f" ({campaign})" if multi else ""))
                    if metrics:
                        counted.setdefault(cfg, []).append(metrics)

                    configs.setdefault(cfg, []).append(run)
                elif os.path.basename(path) == "index.cost.tsv":
                    cost = load(path)
                    if "sshash.index_bytes" in cost:    # Its last line; a failed build leaves it out.
                        costs[(wl, campaign)] = cost

    if not configs:
        sys.exit("No runs found.")

    for campaign, counts in sorted(versions.items()):
        if len(counts) > 1:
            warnings.append(f"{campaign}: runs of {len(counts)} versions: "
                            + ", ".join(f"{c} of {v}" for v, c in sorted(counts.items())))

    for wl in sorted({cfg[0] for cfg in configs}):
        runs = [(cfg, run) for cfg, R in configs.items() if cfg[0] == wl for run in R]
        keys = {"table.size"} | {key for _, run in runs for key in run if re.fullmatch(r"query\..+\.hits", key)}
        for key in sorted(keys):
            vals, known_vals = {}, {}
            for cfg, run in runs:
                if key == "table.size" or key in run:   # Runs may query only some inputs, as `QUERIES` selects.
                    (known_vals if known(cfg) else vals).setdefault(run.get(key, "missing"), set()).add(name(cfg))
            if len(vals) > 1:
                errors.append(f"{wl}: {key} differs: " + "; ".join(f"{v} in {where(c)}" for v, c in vals.items()))
            for v, c in known_vals.items():
                if v not in vals and (vals or len(known_vals) > 1):    # Alone, defective tables can only disagree.
                    defects.append(f"{wl}: {key} is {v} in {where(c)}")

        for cfg, run in runs:
            if "build.ops" in run and run["build.ops"] != run["index.kmers"]:
                errors.append(f"{wl}: {name(cfg)} built {run['build.ops']} of {run['index.kmers']} k-mers")

    queries = []    # The query inputs, in order of first appearance.
    for run in (run for R in configs.values() for run in R):
        for key in run:
            m = re.fullmatch(r"query\.(.+)\.mops", key)
            if m and m[1] not in queries:
                queries.append(m[1])

    query_cols = [f"query.{q}.{s}" for q in queries for s in QUERY_MEDIANS]
    with open(os.path.join(out_dirs[0], "summary.tsv"), "w") as f:
        print("\t".join(CONFIG + ["reps", "version", "table.size"] + MEDIANS + query_cols), file=f)
        for cfg, R in sorted(configs.items(), key=lambda c: (c[0][6],) + c[0][:6]):
            version = ",".join(sorted({run["version"] for run in R if "version" in run}))
            cells = [cfg[6]] + [str(c) for c in cfg[:6]] + [str(len(R)), version, R[0].get("table.size", "")]
            print("\t".join(cells + [median(R, key) for key in MEDIANS + query_cols]), file=f)

    with open(os.path.join(out_dirs[0], "build-cost.tsv"), "w") as f:
        print("campaign\tworkload\tstructure\tsizing\tthreads\tseconds\tmax_rss_bytes\tstages", file=f)
        for (wl, campaign), cost in sorted(costs.items(), key=lambda c: (c[0][1], c[0][0])):
            gg, ss = float(cost["ggcat.seconds"]), float(cost["sshash.seconds"])
            rss = max(int(cost["ggcat.max_rss_bytes"]), int(cost["sshash.max_rss_bytes"]))
            print(f"{campaign}\t{wl}\tsshash\tna\t{cost['threads']}\t{gg + ss:.6g}\t{rss}"
                  f"\tggcat {gg:.6g} + sshash {ss:.6g}", file=f)
            for cfg, R in sorted(configs.items()):
                if cfg[0] == wl and cfg[1] != "sshash" and cfg[2] == "set" and cfg[5] == int(cost["threads"]):
                    ld, it, bd = (float(median(R, f"{p}.seconds")) for p in ("load", "init", "build"))
                    rss = median(R, "memory.build_max_rss_bytes")
                    print(f"{cfg[6]}\t{wl}\t{cfg[1]}.{cfg[3]}\t{cfg[4]}\t{cfg[5]}\t{ld + it + bd:.6g}\t{rss}"
                          f"\tload {ld:.6g} + init {it:.6g} + build {bd:.6g}", file=f)

    with open(os.path.join(out_dirs[0], "counters.tsv"), "w") as f:
        print("\t".join(CONFIG + ["phase", "metric", "value", "runs"]), file=f)
        for cfg, M in sorted(counted.items(), key=lambda c: (c[0][6],) + c[0][:6]):
            vals = {}   # Per phase and metric, the values of the runs that counted it.
            for metrics in M:
                for phase, phase_vals in metrics.items():
                    for metric, val in phase_vals.items():
                        vals.setdefault((phase, metric), []).append(val)
            for (phase, metric), V in vals.items():
                cells = [cfg[6]] + [str(c) for c in cfg[:6]] + [phase, metric, f"{statistics.median(V):.6g}"]
                print("\t".join(cells + [str(len(V))]), file=f)

    print(f"Collated {sum(len(R) for R in configs.values())} runs of {len(configs)} configurations"
          f" into {out_dirs[0]}.")
    for defect in defects:
        print(f"Known defect: {defect}")
    for warning in warnings:
        print(f"Warning: {warning}")
    for error in errors:
        print(error, file=sys.stderr)
    if errors:
        sys.exit(f"{len(errors)} disagreements.")


if __name__ == "__main__":
    main()
