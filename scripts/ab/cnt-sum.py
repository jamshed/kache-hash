#!/usr/bin/env python3
# Summarizes the per-phase counter runs of cnt.sh, fetched from the job directory: per operation, the cycles,
# instructions and each raw event, for each variant and phase, averaged over the groups (cycles and instructions).
# Usage: scripts/ab/cnt-sum.py <dir> [variant...]
import sys, os, re, glob, collections

d = sys.argv[1]
want = sys.argv[2:]
rows = collections.defaultdict(lambda: collections.defaultdict(list))  # (variant, phase) -> event -> values per op
for f in sorted(glob.glob(os.path.join(d, "*.g*.tsv"))):
    v = os.path.basename(f).split(".g")[0]
    if want and v not in want:
        continue
    ops = {}
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        if len(p) != 2:
            continue
        key, val = p
        m = re.match(r"^(\w+)\.ops$", key)
        if m:
            ops[m.group(1)] = float(val)
            continue
        m = re.match(r"^perf\.(\w+)\.(.+)$", key)
        if m and m.group(1) in ops and m.group(2) not in ("task_clock_ns", "page_faults", "context_switches", "cpu_migrations", "running_pct"):
            rows[(v, m.group(1))][m.group(2)].append(float(val) / ops[m.group(1)])

for (v, ph), ev in sorted(rows.items()):
    cyc = sum(ev["cycles"]) / len(ev["cycles"])
    ins = sum(ev["instructions"]) / len(ev["instructions"])
    print(f"== {v} {ph}: {cyc:.1f} cycles/op, {ins:.1f} instr/op, IPC {ins / cyc:.2f}")
    for e, vals in ev.items():
        if e in ("cycles", "instructions"):
            continue
        m = sum(vals) / len(vals)
        share = f"  ({100 * m / cyc:.1f}% of cycles)" if ("stall" in e or "cycles" in e or "walk_active" in e or "pending_cycles" in e or "bound_on" in e or "ports_util" in e or "fb_full" in e) else ""
        print(f"   {e:48s} {m:10.3f}{share}")
