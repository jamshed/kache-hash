# Install rust nightly from release date

`rustup install nightly-2024-09-30`


# Build using nighlty 

`K=31 PREFIX_BITS=24 cargo +nightly-2024-09-30 b --release`

change K and PREFIX_BITS if you want & recompile

# Usage


```
$ ./target/release/bench_cbl -h

A CLI tool with build and query commands

Usage: bench_cbl <COMMAND>

Commands:
  build  Build something from input to output
  query  Query an index with a query file
  help   Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```
