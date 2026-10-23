#!/bin/bash

# Fetches the measurement datasets into a directory: the full set, or a small one for local runs.

data_dir="${1}"
set_name="${2:-"full"}"

if [ $# -lt 1 ] || { [ "${set_name}" != "full" ] && [ "${set_name}" != "small" ]; }; then
    echo "Usage: $0 <data-dir> <(optional) small>" >&2
    exit 1
fi

ncbi="https://ftp.ncbi.nlm.nih.gov/genomes/all"
ucsc="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
ena="https://ftp.sra.ebi.ac.uk/vol1/fastq"


# Downloads the gzipped file at the URL $1 into $2, resuming a partial download.
fetch()
{
    local url="${1}" file="${data_dir}/${2}"

    if [ -s "${file}" ] && gzip -t "${file}" 2> /dev/null; then
        echo "Found ${file}"
        return
    fi

    echo "Fetching ${url}"
    if ! curl -fL --retry 3 -C - -o "${file}" "${url}" || ! gzip -t "${file}"; then
        echo "Failed fetching ${file}" >&2
        rm -f "${file}"
        exit 1
    fi
}


# Downloads the first $3 reads of the gzipped FASTQ at the URL $1 into $2.
fetch_reads()
{
    local url="${1}" file="${data_dir}/${2}" lines=$(( 4 * ${3} ))

    if [ -s "${file}" ] && gzip -t "${file}" 2> /dev/null; then
        echo "Found ${file}"
        return
    fi

    echo "Fetching the first ${3} reads of ${url}"
    curl -fsL "${url}" | gzip -dc 2> /dev/null | head -n "${lines}" | gzip > "${file}"  # The download stops once `head` is done.
    if (( $(gzip -dc "${file}" | wc -l) != lines )); then
        echo "Failed fetching ${file}" >&2
        rm -f "${file}"
        exit 1
    fi
}


mkdir -p "${data_dir}"

fetch "${ncbi}/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" ecoli.fa.gz

if [ "${set_name}" == "small" ]; then
    fetch "${ucsc}/chr21.fa.gz" chr21.fa.gz
    fetch "${ucsc}/chr1.fa.gz" chr1.fa.gz
    fetch_reads "${ena}/SRR583/004/SRR5833294/SRR5833294.fastq.gz" SRR5833294.1M.fastq.gz 1000000
    fetch_reads "${ena}/SRR590/005/SRR5901135/SRR5901135_1.fastq.gz" SRR5901135_1.500K.fastq.gz 500000
else
    fetch "${ncbi}/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" GRCh38.fa.gz
    fetch "${ena}/SRR583/004/SRR5833294/SRR5833294.fastq.gz" SRR5833294.fastq.gz
    fetch "${ena}/SRR590/005/SRR5901135/SRR5901135_1.fastq.gz" SRR5901135_1.fastq.gz
    fetch "${ena}/SRR590/005/SRR5901135/SRR5901135_2.fastq.gz" SRR5901135_2.fastq.gz
fi

ls -l "${data_dir}"
