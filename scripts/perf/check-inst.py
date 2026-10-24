#!/usr/bin/env python3
"""Checks the invariants of the operation statistics in the outputs of instrumented driver runs (`INSTRUMENT=ON`): that
every operation concludes once, with outcomes matching the table and the query hits; that runs, lines, and locks add up;
and that the table's occupancy matches its size. Prints the failures per output, and exits nonzero if any."""

import sys

LEVELS = ["primary", "secondary", "overflow"]
OUTCOMES = ["found", "placed", "absent"]
LINES = ["meta_lines", "slot_read_lines", "slot_write_lines", "overflow_lines"]


def load(path):
    vals = {}
    with open(path) as f:
        for line in f:
            key, _, val = line.rstrip("\n").partition("\t")
            vals[key] = val
    return vals


def check(vals):
    """Returns the failed invariants of the driver output `vals`."""
    fails = []
    num = lambda key: int(vals[key])

    def need(cond, what):
        if not cond:
            fails.append(what)

    concurrent, table_lock = vals["concurrent"] == "1", "no_table_lock" not in vals
    phases = [key[:-4] for key in vals if key.endswith(".ops") and not key.startswith("inst.")]
    for phase in phases:
        prefix, op = f"inst.{phase}.", ("update" if phase == "build" else "query")
        ops = num(prefix + op + ".ops")
        c = {(l, o): num(f"{prefix}{op}.{l}.{o}") for l in LEVELS for o in OUTCOMES}
        found = sum(c[(l, "found")] for l in LEVELS)
        need(ops == num(phase + ".ops"), f"{phase}: {ops} ops counted, {vals[phase + '.ops']} run")
        need(sum(c.values()) == ops, f"{phase}: {sum(c.values())} ops concluded of {ops}")
        need(num(prefix + ("query" if op == "update" else "update") + ".ops") == 0, f"{phase}: ops of the other type")
        need(num(prefix + op + ".key_compares") - num(prefix + op + ".checksum_false_matches") == c[("primary", "found")] + c[("secondary", "found")],
             f"{phase}: equal key compares differ from the main-table finds")
        need(sum(c[("overflow", o)] for o in OUTCOMES) == num(prefix + op + ".overflow_visits"), f"{phase}: overflow visits differ from overflow ops")
        need(num(prefix + op + ".meta_lines") >= ops, f"{phase}: fewer metadata lines than ops")

        need(num(prefix + "runs.len_sum") == ops, f"{phase}: run lengths sum to {vals[prefix + 'runs.len_sum']}, not the ops")
        count = {int(key[len(prefix + "runs.len."):].rstrip("+")): int(val) for key, val in vals.items() if key.startswith(prefix + "runs.len.")}
        need(sum(count.values()) == num(prefix + "runs.count"), f"{phase}: run histogram differs from the run count")
        if max(count, default=0) < 256:     # The last bin also holds the longer runs.
            need(sum(r * n for r, n in count.items()) == ops, f"{phase}: run histogram doesn't sum to the ops")
            need(sum(r * r * n for r, n in count.items()) == num(prefix + "runs.len_sq_sum"), f"{phase}: run histogram doesn't sum to the squared lengths")
        lines = sum(int(val) for key, val in vals.items() if key.startswith(prefix + "runs.lines."))
        need(lines == sum(num(prefix + "runs." + kind) for kind in LINES), f"{phase}: lines by run length don't sum to the runs' lines")

        locks = {key: num(prefix + key) for key in ("bucket_lock.tries", "bucket_lock.fails", "reader_lock.acquisitions", "writer_lock.acquisitions", "resize.count")}
        if op == "update":
            placed = sum(c[(l, "placed")] for l in LEVELS)
            deep = sum(c[(l, o)] for l in ("secondary", "overflow") for o in OUTCOMES)
            acquired = locks["bucket_lock.tries"] - locks["bucket_lock.fails"]
            need(placed == num("table.size"), f"{phase}: {placed} keys placed, table size {vals['table.size']}")
            need(sum(c[(l, "absent")] for l in LEVELS) == 0, f"{phase}: updates without an outcome")
            need(locks["reader_lock.acquisitions"] == (ops if concurrent and table_lock else 0), f"{phase}: {locks['reader_lock.acquisitions']} reader locks")
            if concurrent:
                need(ops + deep <= acquired <= ops + 2 * deep, f"{phase}: {acquired} bucket locks, outside [{ops + deep}, {ops + 2 * deep}]")
            else:
                need(locks["bucket_lock.tries"] == 0, f"{phase}: bucket locks in a sequential table")
            need(num(prefix + "update.overflow_locks") >= c[("overflow", "placed")], f"{phase}: fewer overflow locks than placements")
            need(num("table.capacity") == num("table.init_capacity") << locks["resize.count"], f"{phase}: capacity isn't the initial one doubled per resize")
            need(locks["writer_lock.acquisitions"] >= locks["resize.count"], f"{phase}: fewer writer locks than resizes")
            need((locks["resize.count"] > 0) == (sum(num(prefix + "resize.placed." + l) for l in LEVELS) > 0), f"{phase}: resize placements without resizes")
        else:
            need(found == num(phase + ".hits"), f"{phase}: {found} keys found, {vals[phase + '.hits']} hits")
            need(sum(c[(l, "placed")] for l in LEVELS) == 0, f"{phase}: queries placing keys")
            need(all(v == 0 for v in locks.values()), f"{phase}: queries taking table or bucket locks")
            need(num(prefix + "query.overflow_locks") == c[("overflow", "found")], f"{phase}: overflow locks differ from overflow finds")

    occupancy = {int(key.split(".")[-1]): int(val) for key, val in vals.items() if key.startswith("table.occupancy.")}
    if occupancy:
        need(sum(occupancy.values()) * max(occupancy) == num("table.capacity"), "table: bucket count isn't the capacity over the bucket size")
        need(sum(s * n for s, n in occupancy.items()) == num("table.size") - num("table.overflow_size"), "table: bucket occupancy differs from the main-table size")
    return fails


def main():
    if len(sys.argv) < 2:
        sys.exit(f"Usage: {sys.argv[0]} <driver output TSV>...")

    failed = False
    for path in sys.argv[1:]:
        vals = load(path)
        if vals.get("instrument") != "1":
            print(f"{path}: not an instrumented run")
            failed = True
            continue

        fails = check(vals)
        print(f"{path}: {'OK' if not fails else 'FAILED'}")
        for what in fails:
            print(f"    {what}")
        failed = failed or bool(fails)

    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
