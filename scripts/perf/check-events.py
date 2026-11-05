#!/usr/bin/env python3
"""Checks the raw codes of an events file against `perf list --details` on the current machine."""

import re
import subprocess
import sys


def encoding(name):
    """Returns the raw config of the event `name` per `perf list`, or `None` if perf doesn't know it."""
    out = subprocess.run(["perf", "list", "--details", name], capture_output=True, text=True).stdout
    m = re.search(r"cpu/([^/]*)/", out)
    if not m:
        return None

    fields = dict((kv.split("=") + ["1"])[:2] for kv in m.group(1).split(","))
    val = lambda key: int(fields.get(key, "0"), 0)
    return val("event") | val("umask") << 8 | val("edge") << 18 | val("inv") << 23 | val("cmask") << 24


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <events-file>")

    bad = 0
    for line in open(sys.argv[1]):
        if not line.strip() or line.lstrip().startswith("#"):
            continue

        group, spec = line.split()
        for item in spec.split(","):
            name, code = item.split("=")
            want = encoding(name)
            ok = (want == int(code[1:], 16))
            bad += not ok
            print(f"{group}\t{name}\t{code}\t{'ok' if ok else 'perf: ' + ('unknown' if want is None else f'r{want:x}')}")

    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
