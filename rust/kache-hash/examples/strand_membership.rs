//! Tells, for each read of a FASTA or FASTQ file, how many of its 31-mers occur in a reference's forward strand, and how
//! many in its reverse strand, through a map whose values flag the orientations each k-mer occurs in in the reference.
//! Prints a line per read: its index, its number of k-mer positions, and its k-mers in each strand.
//!
//! Usage: cargo run --release --example strand_membership -- <reference.fa> <reads.fa|reads.fq>

use kache_hash::{Builder, Op, Options};
use std::{env, error::Error, fs, io};

const K: usize = 31;

/// A k-mer's flags: it occurs in the reference as is, i.e. in its canonical form, or as its reverse complement.
const AS_IS: u8 = 1;
const REVCOMP: u8 = 2;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    let (Some(reference), Some(reads)) = (args.get(1), args.get(2)) else {
        return Err("usage: strand_membership <reference.fa> <reads.fa|reads.fq>".into());
    };
    let (reference, reads) = (read_seqs(reference)?, read_seqs(reads)?);

    let positions: u64 = reference.iter().map(|s| (s.len() + 1).saturating_sub(K) as u64).sum();
    let builder = Builder::<u8>::new(&Options { expected_keys: positions.max(1), ..Options::new(K) })?;
    builder.inserter()?.insert_oriented(&reference, Op::Or, AS_IS, REVCOMP)?;
    let flags = builder.freeze();

    // A read's k-mer occurs in the reference's forward strand if it occurs in the reference in its orientation in the
    // read: as is in both, or as the reverse complement in both; and in the reverse strand if in the other one.
    let (mut forward, mut reverse) = (vec![0; reads.len()], vec![0; reads.len()]);
    flags.count_masked(&reads, AS_IS, REVCOMP, &mut forward)?;
    flags.count_masked(&reads, REVCOMP, AS_IS, &mut reverse)?;

    println!("read\tpositions\tforward\treverse");
    for (i, read) in reads.iter().enumerate() {
        println!("{i}\t{}\t{}\t{}", flags.positions(&[read]), forward[i], reverse[i]);
    }
    Ok(())
}

/// Returns the sequences of the uncompressed FASTA or FASTQ file `path`.
fn read_seqs(path: &str) -> io::Result<Vec<Vec<u8>>> {
    let text = fs::read(path)?;
    let lines = text.split(|&c| c == b'\n').map(<[u8]>::trim_ascii_end);
    if text.starts_with(b"@") {
        // FASTQ: the second line of each record of four.
        return Ok(lines.skip(1).step_by(4).map(<[u8]>::to_vec).collect());
    }

    // FASTA: the lines of each record after its header, joined.
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    for line in lines {
        if line.starts_with(b">") {
            seqs.push(Vec::new());
        } else if let Some(s) = seqs.last_mut() {
            s.extend_from_slice(line);
        }
    }
    Ok(seqs)
}
