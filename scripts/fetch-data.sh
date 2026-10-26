#!/bin/bash

# Fetches the measurement and benchmark datasets into a directory: the full set, or a small one for local runs.

data_dir="${1}"
set_name="${2:-"full"}"

if [ $# -lt 1 ] || { [ "${set_name}" != "full" ] && [ "${set_name}" != "small" ]; }; then
    echo "Usage: $0 <data-dir> <(optional) small>" >&2
    exit 1
fi

ncbi="https://ftp.ncbi.nlm.nih.gov/genomes/all"
ucsc="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
ena="https://ftp.sra.ebi.ac.uk/vol1/fastq"
collection="$(dirname "${0}")/ecoli-collection.txt" # NCBI directories of the E. coli collection's genomes.


# Downloads the gzipped file at the URL $1 into $2, resuming a partial download.
fetch()
{
    local url="${1}" file="${data_dir}/${2}" attempt

    if [ -s "${file}" ] && gzip -t "${file}" 2> /dev/null; then
        echo "Found ${file}"
        return
    fi

    echo "Fetching ${url}"
    for attempt in {1..10}; do  # Servers may drop long transfers; each attempt resumes the partial download.
        if curl -fL -C - -o "${file}" "${url}"; then
            if gzip -t "${file}"; then
                return
            fi

            rm -f "${file}"
        fi

        sleep 10
    done

    echo "Failed fetching ${file}" >&2
    exit 1
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


# Downloads the first $1 genomes of the E. coli collection into a subdirectory, and concatenates them into
# `ecoli-$1.fa.gz`: a multi-member gzip, read as one input.
fetch_collection()
{
    local n="${1}" file="${data_dir}/ecoli-${1}.fa.gz" paths path

    if [ -s "${file}" ] && gzip -t "${file}" 2> /dev/null; then
        echo "Found ${file}"
        return
    fi

    paths=$(grep -v '^#' "${collection}" | head -n "${n}")
    mkdir -p "${data_dir}/ecoli-collection"
    for path in ${paths}; do
        fetch "${path}/${path##*/}_genomic.fna.gz" "ecoli-collection/${path##*/}.fna.gz"
    done

    echo "Concatenating ${n} genomes into ${file}"
    for path in ${paths}; do
        cat "${data_dir}/ecoli-collection/${path##*/}.fna.gz"
    done > "${file}.tmp" && mv "${file}.tmp" "${file}"  # A partial concatenation would pass the check above.
}


mkdir -p "${data_dir}"

fetch "${ncbi}/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" ecoli.fa.gz

if [ "${set_name}" == "small" ]; then
    fetch "${ucsc}/chr21.fa.gz" chr21.fa.gz
    fetch "${ucsc}/chr1.fa.gz" chr1.fa.gz
    fetch_reads "${ena}/SRR583/004/SRR5833294/SRR5833294.fastq.gz" SRR5833294.1M.fastq.gz 1000000
    fetch_reads "${ena}/SRR590/005/SRR5901135/SRR5901135_1.fastq.gz" SRR5901135_1.500K.fastq.gz 500000
    fetch_collection 32
else
    fetch "${ncbi}/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" GRCh38.fa.gz
    fetch "${ena}/SRR583/004/SRR5833294/SRR5833294.fastq.gz" SRR5833294.fastq.gz
    fetch "${ena}/SRR590/005/SRR5901135/SRR5901135_1.fastq.gz" SRR5901135_1.fastq.gz
    fetch "${ena}/SRR590/005/SRR5901135/SRR5901135_2.fastq.gz" SRR5901135_2.fastq.gz
    fetch_collection 1024
fi

ls -l "${data_dir}"
