use cbl::CBL;
use needletail::parse_fastx_file;
use std::env::args;
use std::time::{Duration, Instant};

// define the parameters K and T
const K: usize = 25;
type T = u64; // T must be large enough to store $2k + \lg(2k)$ bits
const PREFIX_BITS: usize = 24;

fn main() {
    let args: Vec<String> = args().collect();
    let input_filename = args.get(1).expect("No argument given");

    // create a CBL index with parameters K and T
    let mut cbl = CBL::<K, T, PREFIX_BITS>::new();

    let mut reader = parse_fastx_file(input_filename).unwrap();

    let mut tot_kmer = 0_usize;
    let timer = Instant::now();
    // for each sequence of the FASTA/Q file
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid record");

        let nkmer = if seqrec.seq().len() > K {
            // insert each k-mer of the sequence in the index
            cbl.insert_seq(&seqrec.seq());
            seqrec.seq().len() - K + 1
        } else {
            0
        };
        tot_kmer += nkmer;
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
}
