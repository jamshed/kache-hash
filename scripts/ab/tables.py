#!/usr/bin/env python3
"""Prints a Markdown table of the median Mops/s per configuration, phase and variant of A/B job outputs, with ratios.
The output directories are ab-run.sh's, fetched from the job directories.
Usage: scripts/ab/tables.py <out-dir>[,<out-dir>...] <variant>[,<variant>...] [<a>/<b>[,<a>/<b>...]]"""

import glob
import os
import statistics
import sys


def reps(path):
    """Returns the per-repetition key-value dicts of a driver output."""
    out, cur, head = [], None, {}
    with open(path) as f:
        for line in f:
            key, _, val = line.rstrip("\n").partition("\t")
            if key == "rep":
                cur = dict(head)
                out.append(cur)
            elif cur is None:
                head[key] = val
            else:
                cur[key] = val
    return out or [head]


def main():
    out_dirs, variants = sys.argv[1].split(","), sys.argv[2].split(",")
    ratios = [tuple(x.split("/")) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else []
    data, order = {}, []
    for d in out_dirs:
        for cfg_dir in sorted(glob.glob(os.path.join(d, "*"))):
            if not os.path.isdir(cfg_dir):
                continue
            cfg = os.path.basename(cfg_dir)
            if cfg not in order:
                order.append(cfg)
            for p in sorted(glob.glob(os.path.join(cfg_dir, "*.tsv"))):
                v = os.path.basename(p).rsplit(".r", 1)[0]
                data.setdefault(cfg, {}).setdefault(v, []).extend(reps(p))

    def med(R, key):
        vals = [float(r[key]) for r in R if key in r]
        return statistics.median(vals) if vals else None

    print("| Configuration | Phase | " + " | ".join(variants) + "".join(f" | {a} / {b}" for a, b in ratios) + " |")
    print("|" + "---|" * (2 + len(variants) + len(ratios)))
    for cfg in order:
        for key, name in (("build.mops", "build"), ("query0.mops", "+query"), ("query1.mops", "−query")):
            m = {v: med(data[cfg].get(v, []), key) for v in variants}
            cells = [f"{m[v]:.1f}" if m[v] is not None else "—" for v in variants]
            rs = [f"{m[a] / m[b]:.2f}×" if m.get(a) and m.get(b) else "—" for a, b in ratios]
            print(f"| {cfg if name == 'build' else ''} | {name} | " + " | ".join(cells) + "".join(f" | {r}" for r in rs) + " |")


if __name__ == "__main__":
    main()
