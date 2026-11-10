
# The k-mer counting library: a parallel FASTA and FASTQ reader, and a counter of reads' k-mers within a set, on the
# table. Its interface is plain C++, so that projects linking it compile none of the table's code and keep their own
# flags.

kache_hash_add_library(kache_counter Fastx_Reader.cpp Kmer_Set_Counter.cpp)
