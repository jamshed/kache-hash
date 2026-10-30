#!/bin/bash

# Fetches the benchmark and measurement datasets into a directory: the paper's (`full`), the earlier campaigns'
# (`legacy`), or a small set for local runs and campaign gates (`small`). Read sets are prefixes of runs, with both
# mates, one after the other. HumGut is no longer served, so `full` takes the genomes listed in the file `HUMGUT_LIST`
# from a copy, the directory `HUMGUT`. A new file listed in `data.sha256` must match the sum of its content there.

data_dir="${1}"
set_name="${2:-"full"}"

if [ $# -lt 1 ] || [[ " full legacy small " != *" ${set_name} "* ]]; then
    echo "Usage: $0 <data-dir> <(optional) full | legacy | small>" >&2
    exit 1
fi

ncbi="https://ftp.ncbi.nlm.nih.gov/genomes/all"
ucsc="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
ena="https://ftp.sra.ebi.ac.uk/vol1/fastq"
hpp="https://s3-us-west-2.amazonaws.com/human-pangenomics"
zenodo="https://zenodo.org/api/records"
scripts="$(dirname "${0}")"
collection="${scripts}/ecoli-collection.txt" # NCBI directories of the E. coli collection's genomes.
manifest="${scripts}/data.sha256"   # SHA-256 sums of the full set's files' contents, as fetched for the paper.


# Prints the SHA-256 sum of the content of the gzipped file $1: its bytes vary with the gzip that wrote them.
sha256_of()
{
    gzip -dc "${1}" | { sha256sum 2> /dev/null || shasum -a 256; } | cut -d' ' -f1
}


# Prints the MD5 sum of the file $1.
md5_of()
{
    { md5sum "${1}" 2> /dev/null || md5 -r "${1}"; } | cut -d' ' -f1
}


# Checks the new file $1 of the data directory against its sum in the manifest, if listed there; else prints its sum,
# for the manifest. Exits if they differ.
accept()
{
    local file="${data_dir}/${1}" sum

    [[ "${1}" == */* ]] && return   # The members of a collection.
    sum=$(sha256_of "${file}")
    if grep -q "  ${1}$" "${manifest}" 2> /dev/null; then
        if ! grep -q "^${sum}  ${1}$" "${manifest}"; then
            echo "${file} differs from its sum in ${manifest}" >&2
            exit 1
        fi
    else
        echo "Not in ${manifest}: ${sum}  ${1}"
    fi
}


# Downloads the gzipped file at the URL $1 into $2, resuming a partial download, and checks it against the MD5 sum $3
# that its publisher gives, if any.
fetch()
{
    local url="${1}" file="${data_dir}/${2}" md5="${3:-}" attempt

    if [ -s "${file}" ] && gzip -t "${file}" 2> /dev/null; then
        echo "Found ${file}"
        return
    fi

    echo "Fetching ${url}"
    for attempt in {1..10}; do  # Servers may drop long transfers; each attempt resumes the partial download.
        if curl -fL -C - -o "${file}" "${url}"; then
            if gzip -t "${file}"; then
                if [ -n "${md5}" ] && [ "$(md5_of "${file}")" != "${md5}" ]; then
                    echo "${file} differs from its published MD5 sum" >&2
                    rm -f "${file}"
                    exit 1
                fi

                accept "${2}"
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
    curl -fsL "${url}" | gzip -dc 2> /dev/null | head -n "${lines}" | gzip -n > "${file}"  # The download stops once `head` is done.
    if (( $(gzip -dc "${file}" | wc -l) != lines )); then
        echo "Failed fetching ${file}" >&2
        rm -f "${file}"
        exit 1
    fi

    accept "${2}"
}


# Writes the first $3 reads of the gzipped FASTQ $2 of the data directory into $1.
prefix_reads()
{
    local file="${data_dir}/${1}" lines=$(( 4 * ${3} ))

    if [ -s "${file}" ] && gzip -t "${file}" 2> /dev/null; then
        echo "Found ${file}"
        return
    fi

    echo "Writing the first ${3} reads of ${2} into ${file}"
    gzip -dc "${data_dir}/${2}" | head -n "${lines}" | gzip -n > "${file}.tmp"
    if (( $(gzip -dc "${file}.tmp" | wc -l) != lines )); then
        echo "Failed writing ${file}" >&2
        rm -f "${file}.tmp"
        exit 1
    fi

    mv "${file}.tmp" "${file}"
    accept "${1}"
}


# Concatenates the gzipped files $2..., in order, into $1: a multi-member gzip, read as one input. Paths are in the
# data directory unless absolute.
concat()
{
    local file="${data_dir}/${1}" parts=() part
    shift

    if [ -s "${file}" ] && gzip -t "${file}" 2> /dev/null; then
        echo "Found ${file}"
        return
    fi

    for part in "$@"; do
        [[ "${part}" == /* ]] || part="${data_dir}/${part}"
        if [ ! -s "${part}" ]; then
            echo "Missing ${part} for ${file}" >&2
            exit 1
        fi

        parts+=("${part}")
    done

    echo "Concatenating ${#parts[@]} files into ${file}"   # Through a temporary file: a partial one would be found.
    cat "${parts[@]}" > "${file}.tmp" && mv "${file}.tmp" "${file}"
    accept "$(basename "${file}")"
}


# Downloads the first $1 genomes of the E. coli collection into a subdirectory, and concatenates them into
# `ecoli-$1.fa.gz`.
fetch_collection()
{
    local n="${1}" file="${data_dir}/ecoli-${1}.fa.gz" paths path parts=()

    if [ -s "${file}" ] && gzip -t "${file}" 2> /dev/null; then
        echo "Found ${file}"
        return
    fi

    paths=$(grep -v '^#' "${collection}" | head -n "${n}")
    mkdir -p "${data_dir}/ecoli-collection"
    for path in ${paths}; do
        fetch "${path}/${path##*/}_genomic.fna.gz" "ecoli-collection/${path##*/}.fna.gz"
        parts+=("ecoli-collection/${path##*/}.fna.gz")
    done

    concat "ecoli-${n}.fa.gz" "${parts[@]}"
}


# Downloads the 3,682 E. coli assemblies of Themisto and Fulgor (Zenodo record 6577997), and concatenates them in
# file-name order into `ecoli-3682.fa.gz`, and every 8th of them into `ecoli-3682.every8.fa.gz`.
fetch_ecoli_3682()
{
    local dir="${data_dir}/coli3682_dataset" cpus genomes

    if [ -s "${data_dir}/ecoli-3682.fa.gz" ] && [ -s "${data_dir}/ecoli-3682.every8.fa.gz" ]; then
        echo "Found ${data_dir}/ecoli-3682.fa.gz and ${data_dir}/ecoli-3682.every8.fa.gz"
        return
    fi

    fetch "${zenodo}/6577997/files/coli3682_dataset.tar.gz/content" coli3682_dataset.tar.gz \
        244cde0be80717d2a7a0af0c47eb9f86
    if [ ! -e "${dir}/.extracted" ]; then
        echo "Extracting ${data_dir}/coli3682_dataset.tar.gz"
        tar -xzf "${data_dir}/coli3682_dataset.tar.gz" -C "${data_dir}" && touch "${dir}/.extracted" || exit 1
    fi

    cpus=$(nproc 2> /dev/null || sysctl -n hw.ncpu)    # The collection's members get gzipped, in parallel.
    find "${dir}" -name "*.fna" -print0 | xargs -0 -r -n 16 -P "${cpus}" gzip -n -f
    genomes=$(cd "${dir}" && ls | grep '\.fna\.gz$' | LC_ALL=C sort)
    if (( $(echo "${genomes}" | wc -l) != 3682 )); then
        echo "Found $(echo "${genomes}" | wc -l) genomes in ${dir}, not 3,682" >&2
        exit 1
    fi

    concat ecoli-3682.fa.gz $(echo "${genomes}" | sed 's|^|coli3682_dataset/|')
    concat ecoli-3682.every8.fa.gz $(echo "${genomes}" | awk 'NR % 8 == 1' | sed 's|^|coli3682_dataset/|')
}


# Concatenates the HumGut genomes listed in the file `HUMGUT_LIST`, in order, from a copy of HumGut's per-genome files
# in the directory `HUMGUT`, into `humgut-1in20.fa.gz`, and every 8th of them into `humgut-1in20.every8.fa.gz`.
fetch_humgut()
{
    local names

    if [ -s "${data_dir}/humgut-1in20.fa.gz" ] && [ -s "${data_dir}/humgut-1in20.every8.fa.gz" ]; then
        echo "Found ${data_dir}/humgut-1in20.fa.gz and ${data_dir}/humgut-1in20.every8.fa.gz"
        return
    fi

    if [ -z "${HUMGUT:-}" ] || [ ! -s "${HUMGUT_LIST:-}" ]; then
        echo "Set HUMGUT to a directory of HumGut's genomes, which is no longer served, and HUMGUT_LIST to a list" >&2
        exit 1
    fi

    names=$(grep -v '^#' "${HUMGUT_LIST}")
    concat humgut-1in20.fa.gz $(echo "${names}" | sed "s|^|${HUMGUT}/|")
    concat humgut-1in20.every8.fa.gz $(echo "${names}" | awk 'NR % 8 == 1' | sed "s|^|${HUMGUT}/|")
}


mkdir -p "${data_dir}"

fetch "${ncbi}/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" ecoli.fa.gz

if [ "${set_name}" == "small" ]; then
    fetch "${ucsc}/chr21.fa.gz" chr21.fa.gz
    fetch "${ucsc}/chr1.fa.gz" chr1.fa.gz
    fetch_reads "${ena}/SRR583/004/SRR5833294/SRR5833294.fastq.gz" SRR5833294.1M.fastq.gz 1000000
    fetch_reads "${ena}/SRR590/005/SRR5901135/SRR5901135_1.fastq.gz" SRR5901135_1.500K.fastq.gz 500000
    fetch_collection 32
elif [ "${set_name}" == "legacy" ]; then
    fetch "${ncbi}/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" GRCh38.fa.gz
    fetch "${ena}/SRR583/004/SRR5833294/SRR5833294.fastq.gz" SRR5833294.fastq.gz
    fetch "${ena}/SRR590/005/SRR5901135/SRR5901135_1.fastq.gz" SRR5901135_1.fastq.gz
    fetch "${ena}/SRR590/005/SRR5901135/SRR5901135_2.fastq.gz" SRR5901135_2.fastq.gz
    fetch_collection 1024
else
    # The reference: T2T-CHM13v2.0, the analysis set (soft-masked, with chrY and chrM).
    fetch "${hpp}/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz" CHM13.fa.gz 06caded0055ec647a69acb708e13beff

    # Short reads: HG002, NovaSeq X 2x150 PCR-free, the first 20M pairs, about 2x; and their first 4M pairs, a query.
    fetch_reads "${ena}/SRR373/039/SRR37356339/SRR37356339_1.fastq.gz" SRR37356339_1.20M.fastq.gz 20000000
    fetch_reads "${ena}/SRR373/039/SRR37356339/SRR37356339_2.fastq.gz" SRR37356339_2.20M.fastq.gz 20000000
    concat SRR37356339.20M.fastq.gz SRR37356339_1.20M.fastq.gz SRR37356339_2.20M.fastq.gz
    prefix_reads SRR37356339_1.4M.fastq.gz SRR37356339_1.20M.fastq.gz 4000000
    prefix_reads SRR37356339_2.4M.fastq.gz SRR37356339_2.20M.fastq.gz 4000000
    concat SRR37356339.4M.fastq.gz SRR37356339_1.4M.fastq.gz SRR37356339_2.4M.fastq.gz

    # Long reads: HG002, PacBio HiFi of a 15-kb Sequel II library, the first 360K reads, about 1.5x; and their first
    # 45K reads, a query.
    fetch_reads "${ena}/SRR103/046/SRR10382246/SRR10382246_subreads.fastq.gz" SRR10382246.360K.fastq.gz 360000
    prefix_reads SRR10382246.45K.fastq.gz SRR10382246.360K.fastq.gz 45000

    # The pangenome, and a query of its every 8th genome.
    fetch_ecoli_3682

    # The metagenome, and a query of its every 8th genome.
    fetch_humgut

    # Further queries: E. coli reads, MiSeq paired-end of about 200 bases, both mates; and a human stool metagenome
    # (MetaHIT), HiSeq 2000 2x90, both mates.
    fetch "${ena}/SRR590/005/SRR5901135/SRR5901135_1.fastq.gz" SRR5901135_1.fastq.gz 45da7d3f77b4ef01cbcd90469f3f8312
    fetch "${ena}/SRR590/005/SRR5901135/SRR5901135_2.fastq.gz" SRR5901135_2.fastq.gz 27ff8143298bc16659bb865e889430ff
    concat SRR5901135.fastq.gz SRR5901135_1.fastq.gz SRR5901135_2.fastq.gz
    fetch "${ena}/ERR321/ERR321482/ERR321482_1.fastq.gz" ERR321482_1.fastq.gz 9a356eb3f58722d049115b95d48402fb
    fetch "${ena}/ERR321/ERR321482/ERR321482_2.fastq.gz" ERR321482_2.fastq.gz 767d684f9f7ebc89b036b665392e4f31
    concat ERR321482.fastq.gz ERR321482_1.fastq.gz ERR321482_2.fastq.gz
fi

ls -l "${data_dir}"
