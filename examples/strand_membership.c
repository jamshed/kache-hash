
/*
 * Tells, for each read of a FASTA or FASTQ file, how many of its 31-mers occur
 * in a reference's forward strand, and how many in its reverse strand, through
 * kache-hash's C interface, as the C++ and Rust examples `strand_membership`
 * do: a map of 8-bit masks of the orientations each k-mer occurs in in the
 * reference, and two masked counts of the reads' k-mers. Prints a line per
 * read: its index, its number of k-mer positions, and its k-mers in each
 * strand. One thread inserts the reference here; each further one would
 * register a user of its own.
 * Usage: strand_membership_c <reference.fa> <reads.fa|reads.fq>
 */

#include <kache_hash/kache_capi.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


enum { K = 31 };

/* A k-mer's flags: it occurs in the reference as is, i.e. in its canonical
 * form, or as its reverse complement. */
enum { AS_IS = 1, REVCOMP = 2 };


/* The sequences of a file: `n` of them at `seqs`, of lengths `lens`, within the
 * file's text at `text`. */
typedef struct
{
    char* text;
    char** seqs;
    size_t* lens;
    size_t n;
} seq_set;


/* Reads the sequences of the uncompressed FASTA or FASTQ file at `path` into
 * `s`. A FASTA record's lines are joined in place. Returns 0 on success. */
static int read_seqs(const char* const path, seq_set* const s)
{
    FILE* const f = fopen(path, "rb");
    if(f == NULL)
        return -1;

    fseek(f, 0, SEEK_END);
    const long size = ftell(f);
    fseek(f, 0, SEEK_SET);
    s->text = malloc((size_t)size + 1);
    if(s->text == NULL || fread(s->text, 1, (size_t)size, f) != (size_t)size)
    {
        fclose(f);
        return -1;
    }
    fclose(f);
    s->text[size] = '\n';

    /* At most a sequence per line. */
    size_t lines = 0;
    for(long i = 0; i <= size; ++i)
        lines += s->text[i] == '\n';
    s->seqs = malloc(lines * sizeof(char*));
    s->lens = malloc(lines * sizeof(size_t));
    if(s->seqs == NULL || s->lens == NULL)
        return -1;

    const int fastq = size > 0 && s->text[0] == '@';
    s->n = 0;
    char* out = NULL;   /* The end of the current FASTA sequence, joined so far. */
    size_t line = 0;
    for(char* p = s->text; p < s->text + size; ++line)
    {
        char* const eol = memchr(p, '\n', (size_t)(s->text + size - p) + 1);
        size_t len = (size_t)(eol - p);
        if(len > 0 && p[len - 1] == '\r')
            len--;

        if(fastq)
        {
            if(line % 4 == 1)
                s->seqs[s->n] = p, s->lens[s->n++] = len;
        }
        else if(len > 0 && p[0] == '>')
            out = NULL;
        else if(out == NULL)
            s->seqs[s->n] = p, s->lens[s->n++] = len, out = p + len;
        else
            memmove(out, p, len), out += len, s->lens[s->n - 1] += len;

        p = eol + 1;
    }

    return 0;
}


/* Exits with the description of the status `code` of the call `what`, unless
 * it is `KACHE_OK`. */
static void check(const int code, const char* const what)
{
    if(code != KACHE_OK)
    {
        fprintf(stderr, "%s: %s\n", what, kache_strerror(code));
        exit(1);
    }
}


int main(const int argc, const char* const argv[])
{
    if(argc < 3)
    {
        fprintf(stderr, "usage: strand_membership_c <reference.fa> <reads.fa|reads.fq>\n");
        return 1;
    }

    seq_set ref, reads;
    if(read_seqs(argv[1], &ref) != 0 || read_seqs(argv[2], &reads) != 0)
    {
        fprintf(stderr, "cannot read the inputs\n");
        return 1;
    }

    /* A map of 8-bit flags, sized for the reference's k-mer positions, which
     * bound its distinct k-mers; it grows beyond as needed. */
    kache_options opt = kache_default_options(K, 1);
    opt.expected_keys = 1;
    for(size_t i = 0; i < ref.n; ++i)
        opt.expected_keys += ref.lens[i] >= K ? ref.lens[i] - K + 1 : 0;

    kache_table* flags;
    check(kache_new(&opt, &flags), "kache_new");
    uint32_t user;
    check(kache_register(flags, &user), "kache_register");
    check(kache_insert_oriented(flags, user, (const char* const*)ref.seqs, ref.lens, ref.n, KACHE_OR, AS_IS, REVCOMP, NULL), "kache_insert_oriented");

    /* A read's k-mer occurs in the reference's forward strand if it occurs in
     * the reference in its orientation in the read: as is in both, or as the
     * reverse complement in both; and in the reverse strand if in the other
     * one. */
    uint64_t* const forward = malloc((reads.n + 1) * sizeof(uint64_t));
    uint64_t* const reverse = malloc((reads.n + 1) * sizeof(uint64_t));
    if(forward == NULL || reverse == NULL)
        return 1;
    check(kache_count_masked(flags, (const char* const*)reads.seqs, reads.lens, reads.n, AS_IS, REVCOMP, forward), "kache_count_masked");
    check(kache_count_masked(flags, (const char* const*)reads.seqs, reads.lens, reads.n, REVCOMP, AS_IS, reverse), "kache_count_masked");

    printf("read\tpositions\tforward\treverse\n");
    for(size_t i = 0; i < reads.n; ++i)
        printf("%zu\t%zu\t%llu\t%llu\n", i, reads.lens[i] >= K ? reads.lens[i] - K + 1 : 0,
               (unsigned long long)forward[i], (unsigned long long)reverse[i]);

    kache_free(flags);
    free(forward), free(reverse);
    free(ref.text), free(ref.seqs), free(ref.lens);
    free(reads.text), free(reads.seqs), free(reads.lens);
    return 0;
}
