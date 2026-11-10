//! Counts the 31-mers of a FASTA or FASTQ file from a number of threads, and prints the number of distinct k-mers, of
//! their occurrences, and of the k-mers by their counts, up to 10.
//!
//! Usage: cargo run --release --example kmer_counts -- <seqs.fa|seqs.fq> [threads]

use kache_hash::{Builder, Op, Options, Reduction};
use std::{env, error::Error, fs, io, thread};

const K: usize = 31;

/// The counts that the histogram tells apart: 1 to `MAX_COUNT - 1`, and `MAX_COUNT` or more.
const MAX_COUNT: usize = 10;

/// The k-mer positions of a piece of a sequence that a thread inserts at once.
const PIECE: usize = 1 << 14;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    let path = args.get(1).ok_or("usage: kmer_counts <seqs.fa|seqs.fq> [threads]")?;
    let threads: usize = args.get(2).map_or(Ok(4), |t| t.parse())?;
    if threads == 0 {
        return Err("threads: at least 1".into());
    }
    let seqs = read_seqs(path)?;

    // The table grows beyond its expected keys as needed; the k-mer positions bound them, and so spare it the resizes.
    let positions: u64 = seqs.iter().map(|s| (s.len() + 1).saturating_sub(K) as u64).sum();
    let builder = Builder::<u32>::new(&Options { expected_keys: positions.max(1), ..Options::new(K) })?;

    // An inserter per thread, each inserting a part of the sequences' pieces, so that the threads share even one genome.
    let pieces = pieces(&seqs);
    thread::scope(|s| {
        let builder = &builder;
        let parts = pieces.chunks(pieces.len().div_ceil(threads).max(1));
        let handles: Vec<_> = parts.map(|part| s.spawn(move || builder.inserter()?.insert(part, Op::Add, 1))).collect();
        handles.into_iter().try_for_each(|h| h.join().unwrap().map(|_| ()))
    })?;
    let counts = builder.freeze();

    // The histogram, from a share of the table per thread.
    let mut hist = [0u64; MAX_COUNT + 1];
    thread::scope(|s| {
        let counts = &counts;
        let handles: Vec<_> = (0..threads)
            .map(|shard| {
                s.spawn(move || {
                    let mut h = [0u64; MAX_COUNT + 1];
                    counts
                        .for_each_chunk(threads, shard, |_, values| {
                            values.iter().for_each(|&c| h[(c as usize).min(MAX_COUNT)] += 1);
                            true
                        })
                        .map(|()| h)
                })
            })
            .collect();
        for handle in handles {
            hist.iter_mut().zip(handle.join().unwrap()?).for_each(|(a, b)| *a += b);
        }
        Ok::<_, kache_hash::Error>(())
    })?;

    println!("distinct {K}-mers\t{}", counts.len());
    println!("occurrences\t{}", counts.reduce(Reduction::Sum, u32::MAX, threads)?);
    for (c, n) in hist.iter().enumerate().skip(1) {
        println!("count {c}{}\t{n}", if c == MAX_COUNT { "+" } else { "" });
    }
    Ok(())
}

/// Returns the pieces of the sequences `seqs` of at most `PIECE` k-mer positions each, overlapping by `K - 1` bases, so
/// that their k-mer positions are the sequences'.
fn pieces(seqs: &[Vec<u8>]) -> Vec<&[u8]> {
    let mut pieces = Vec::new();
    for s in seqs {
        let mut i = 0;
        while i + K <= s.len() {
            pieces.push(&s[i..(i + PIECE + K - 1).min(s.len())]);
            i += PIECE;
        }
    }
    pieces
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
