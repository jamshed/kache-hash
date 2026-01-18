#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

// define the parameters K and T
const K: usize = 31;
type T = u128; // T must be large enough to store $2k + \lg(2k)$ bits
const PREFIX_BITS: usize = 24;

use bincode::{DefaultOptions, Options};
use cbl::cbl::kmer_bits;
use cbl::{kmer::Kmer, CBL};
use clap::{Parser, Subcommand};
use needletail::parse_fastx_file;
use serde::{de::DeserializeOwned, Serialize};
use std::env::args;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

pub struct Timer {
    total: Duration,
    start: Option<Instant>,
}

impl Timer {
    pub fn new() -> Self {
        Self {
            total: Duration::ZERO,
            start: None,
        }
    }

    pub fn start(&mut self) {
        self.start = Some(Instant::now());
    }

    pub fn stop(&mut self) {
        if let Some(start) = self.start.take() {
            self.total += start.elapsed();
        }
    }

    pub fn elapsed(&self) -> Duration {
        self.total
    }

    pub fn elapsed_secs(&self) -> f64 {
        self.total.as_secs_f64()
    }
}

#[derive(Parser)]
#[command(name = "myapp")]
#[command(about = "A CLI tool with build and query commands", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build something from input to output
    Build {
        /// Input file path
        #[arg(short, long)]
        input: PathBuf,

        /// Output file path
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Don't use canonical k-mers (default is canonical)
        #[arg(short, long)]
        non_canonical: bool,
    },
    /// Query an index with a query file
    Query {
        /// Index file path
        #[arg(short, long)]
        index: PathBuf,

        /// Query file path
        #[arg(short, long)]
        query: PathBuf,
    },
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Build {
            input,
            output,
            non_canonical,
        } => {
            let input_filename = input.as_path();

            // create a CBL index with parameters K and T
            let mut cbl = if *non_canonical {
                CBL::<K, T, PREFIX_BITS>::new()
            } else {
                CBL::<K, T, PREFIX_BITS>::new_canonical()
            };

            let mut reader = parse_fastx_file(input_filename).unwrap();

            let mut tot_kmer = 0_usize;
            let mut timer = Timer::new();
            // for each sequence of the FASTA/Q file
            while let Some(record) = reader.next() {
                let seqrec = record.expect("Invalid record");

                let nkmer = if seqrec.seq().len() > K {
                    // insert each k-mer of the sequence in the index
                    timer.start();
                    cbl.insert_seq(&seqrec.seq());
                    timer.stop();
                    seqrec.seq().len() - K + 1
                } else {
                    0
                };
                tot_kmer += nkmer;
            }
            if let Some(output_filename) = output {
                write_index(&cbl, output_filename.as_path());
            }

            let dur = timer.elapsed();
            println!(
                "Time taken for cbl (sequential) for insertion ot {tot_kmer} k-mers was {}s.",
                dur.as_secs_f64()
            );
            println!(
                "Throughput: {} M ops/s",
                (tot_kmer as f64 / (dur.as_secs_f64() * 1_000_000_f64))
            );
            Ok(())
        }
        Commands::Query { index, query } => {
            // Your query logic here
            let index_filename = index.as_path();
            let mut cbl: CBL<K, T, PREFIX_BITS> = read_index(index_filename);

            let mut reader = parse_fastx_file(query.as_path()).unwrap();

            let mut nquery = 0_usize;
            let mut npos = 0_usize;
            let mut nneg = 0_usize;
            let mut timer = Timer::new();
            // for each sequence of the FASTA/Q file
            while let Some(record) = reader.next() {
                let seqrec = record.expect("Invalid record");
                if seqrec.seq().len() > K {
                    // insert each k-mer of the sequence in the index
                    timer.start();
                    let contains_vec = cbl.contains_seq(&seqrec.seq());
                    timer.stop();
                    let pos_queries = contains_vec
                        .iter()
                        .filter(|b| **b) // Filter keeps only elements for which the closure returns true
                        .count();
                    let n = contains_vec.len();
                    npos += pos_queries;
                    nneg += n - pos_queries;
                    nquery += n;
                }
            }

            let dur = timer.elapsed();
            println!(
                "Time taken for cbl (sequential) for query of {nquery} k-mers was {}s.",
                dur.as_secs_f64()
            );
            println!(
                "Throughput: {} M ops/s",
                (nquery as f64 / (dur.as_secs_f64() * 1_000_000_f64))
            );
            println!(
                "positive: {}, npos: {}, hit percentage: {}%",
                npos,
                nneg,
                100f64 * (npos as f64) / (nquery as f64)
            );
            Ok(())
        }
    }
}

fn write_index<S: Serialize, P: AsRef<Path> + Copy>(index: &S, path: P) {
    let output = File::create(path)
        .unwrap_or_else(|_| panic!("Failed to open {}", path.as_ref().to_str().unwrap()));
    let mut writer = BufWriter::new(output);
    eprintln!("Writing the index to {}", path.as_ref().to_str().unwrap());
    DefaultOptions::new()
        .with_varint_encoding()
        .reject_trailing_bytes()
        .serialize_into(&mut writer, &index)
        .unwrap();
}

fn read_index<D: DeserializeOwned, P: AsRef<Path> + Copy>(path: P) -> D {
    let index = File::open(path)
        .unwrap_or_else(|_| panic!("Failed to open {}", path.as_ref().to_str().unwrap()));
    let reader = BufReader::new(index);
    eprintln!(
        "Reading the index stored in {}",
        path.as_ref().to_str().unwrap()
    );
    DefaultOptions::new()
        .with_varint_encoding()
        .reject_trailing_bytes()
        .deserialize_from(reader)
        .unwrap()
}
