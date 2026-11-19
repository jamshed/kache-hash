#!/usr/bin/env python3
# Doxygen's input filter for kache-hash's headers, whose documentation is in plain `//` comments: it turns the
# comments before a declaration into `///` ones, a comment after a member on its line into a `///<` one, with the
# comment-only lines that continue it, and a file's first comment, before any declaration and followed by a blank
# line, into the file's own. Comments in function bodies stay plain, so that they stay out of the functions' docs; the
# lines keep their numbers but for the file's comment, which takes one more.
# Usage, from Doxygen's FILTER_PATTERNS: doxygen_filter.py <header>

import re
import sys

# The heads of scopes of declarations, rather than of code: the text before their `{`, with `alignas` and attributes
# taken out, and strings emptied, as `extern "C"`'s.
DECL_HEAD = re.compile(r'^(template\s*<[^{}]*>\s*)?(typedef\s+)?(class|struct|union|enum|namespace)\b[^(){}]*$|^extern\s+""$')
NOISE = re.compile(r'alignas\s*\([^()]*\)|__attribute__\s*\(\([^()]*(\([^()]*\)[^()]*)*\)\)')


def code_and_comment(line):
    """Splits `line` into its code and the column of its `//` comment, if any, outside strings and characters."""
    i, quote = 0, None
    while i < len(line):
        c = line[i]
        if quote:
            if c == '\\':
                i += 1
            elif c == quote:
                quote = None
        elif c in '"\'':
            quote = c
        elif line.startswith('//', i):
            return line[:i], i
        elif line.startswith('/*', i):
            end = line.find('*/', i + 2)
            if end < 0:
                return line[:i], None
            i = end + 1
        i += 1
    return line, None


def strip_literals(code):
    """Returns `code` without the contents of its strings and characters, and its block comments."""
    code = re.sub(r'/\*.*?\*/', ' ', code)
    return re.sub(r'"(\\.|[^"\\])*"|\'(\\.|[^\'\\])*\'', '""', code)


def main(path):
    lines = open(path, encoding='utf-8').read().split('\n')
    out = []
    stack = []  # The open scopes: whether each is one of code.
    head = ''   # The code since the last `;`, `{` or `}`.
    seen_decl = False   # Whether a declaration came before, outside the preprocessor.
    trail_col = None    # The column of the trailing comment that the line before ended with.
    in_macro = False    # Whether the line continues a preprocessor directive.
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.lstrip()
        in_code = any(stack)

        if in_macro or stripped.startswith('#'):
            in_macro = line.endswith('\\')
            code, col = code_and_comment(line)
            if col is not None and not in_code and stripped.startswith('#define') and not line.startswith('///', col):
                line = code + '///<' + line[col + 2:]
            out.append(line)
            trail_col = None
            i += 1
            continue

        if stripped.startswith('//'):
            col = len(line) - len(stripped)
            plain = not stripped.startswith(('///', '//!'))
            if in_code or not plain:
                out.append(line)
            elif trail_col is not None and col == trail_col:
                out.append(line[:col] + '///<' + line[col + 2:])
            elif not seen_decl and not stack:
                # The file's first comment: its own documentation, iff a blank line follows the block.
                j = i
                while j < len(lines) and lines[j].lstrip().startswith('//'):
                    j += 1
                block = [l[:len(l) - len(l.lstrip())] + '///' + l.lstrip()[2:] for l in lines[i:j]]
                if j < len(lines) and lines[j].strip() == '':
                    out.append('/// @file')
                out.extend(block)
                seen_decl = seen_decl or not (j < len(lines) and lines[j].strip() == '')
                i = j
                continue
            else:
                out.append(line[:col] + '///' + line[col + 2:])
            i += 1
            continue

        code, col = code_and_comment(line)
        trail_col = None
        if col is not None and not in_code and not line.startswith('///', col):
            ends = code.rstrip()[-1:]
            if ends and ends not in '{}:':
                line = code + '///<' + line[col + 2:]
                trail_col = col
        out.append(line)

        # Track the scopes that the line's code opens and closes.
        for c in strip_literals(code):
            if c == '{':
                h = NOISE.sub(' ', head).strip()
                h = re.sub(r'\s+', ' ', h)
                stack.append(not DECL_HEAD.match(h))
                head = ''
            elif c == '}':
                if stack:
                    stack.pop()
                head = ''
            elif c == ';':
                head = ''
            else:
                head += c
        head += ' '
        if code.strip():
            seen_decl = True
        i += 1

    sys.stdout.write('\n'.join(out))


if __name__ == '__main__':
    main(sys.argv[1])
