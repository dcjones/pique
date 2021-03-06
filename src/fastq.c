
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "fastq.h"
#include "misc.h"


static void str_init(str_t* str)
{
    str->n = 0;
    str->size = 128;
    str->s = malloc_or_die(str->size);
    str->s[0] = '\0';
}


static void str_free(str_t* str)
{
    free(str->s);
}


/* Reserve space for `size` more characters. */
static void str_reserve_extra(str_t* str, size_t size)
{
    if (str->n + size > str->size) {
        if (str->n + size > 2 * str->size) {
            str->size = str->n + size;
        }
        else str->size *= 2;
        str->s = realloc_or_die(str->s, str->size * sizeof(char));
    }
}


/* Copy n characters from c to the end of str. */
static void str_append(str_t* str, char* c, size_t n)
{
    str_reserve_extra(str, n);
    memcpy(str->s + str->n, c, n);
    str->n += n;
    str->s[str->n] = '\0';
}


seq_t* seq_create()
{
    seq_t* seq = malloc_or_die(sizeof(seq_t));
    str_init(&seq->id1);
    str_init(&seq->seq);
    str_init(&seq->id2);
    str_init(&seq->qual);
    return seq;
}


void seq_free(seq_t* seq)
{
    str_free(&seq->id1);
    str_free(&seq->seq);
    str_free(&seq->id2);
    str_free(&seq->qual);
    free(seq);
}


static const size_t parser_buf_size = 1000000;


struct fastq_t_
{
    FILE* file;
    size_t readlen;
    char* buf;
    char* next;
    bool linestart;
};



fastq_t* fastq_create(FILE* file)
{
    fastq_t* f = malloc_or_die(sizeof(fastq_t));
    f->file = file;
    f->next = f->buf = malloc_or_die(parser_buf_size);
    f->readlen = 0;
    f->linestart = true;
    return f;
}


void fastq_free(fastq_t* f)
{
    free(f->buf);
    free(f);
}



bool fasta_read(fastq_t* f, seq_t* seq)
{
    enum {
        FASTA_STATE_INIT,
        FASTA_STATE_ID,   /* Reading ID */
        FASTA_STATE_SEQ,  /* Reading sequence. */
    } state = FASTA_STATE_INIT;

    seq->id1.n = seq->seq.n = seq->id2.n = seq->qual.n = 0;
    char* end = f->buf + f->readlen;
    do {
        while (f->next < end) {
            if (f->linestart && f->next[0] == '>') {
                if (seq->id1.n > 0) return true;
                state = FASTA_STATE_ID;
                f->linestart = false;
                ++f->next;
                continue;
            }

            char* u = memchr(f->next, '\n', end - f->next);
            if (u == NULL) {
                f->linestart = false;
                u = end;
            }
            else f->linestart = true;

            if (state == FASTA_STATE_ID) {
                str_append(&seq->id1, f->next, u - f->next);
                if (f->linestart) state = FASTA_STATE_SEQ;
            }
            else if (state == FASTA_STATE_SEQ) {
                char* v = memchr(f->next, ' ', end - f->next);
                if (v != NULL && v < u) u = v;
                str_append(&seq->seq, f->next, u - f->next);
            }

            f->next = u + 1;
        }

        /* Try to read more. */
        f->readlen = fread(f->buf, 1, parser_buf_size, f->file);
        f->next = f->buf;
        end = f->buf + f->readlen;
    } while (f->readlen);

    return false;
}


bool fastq_read(fastq_t* f, seq_t* seq)
{
    enum {
        FASTQ_STATE_ID1,  /* Reading ID1. */
        FASTQ_STATE_SEQ,  /* Reading the sequence. */
        FASTQ_STATE_ID2,  /* Reading ID2. */
        FASTQ_STATE_QUAL, /* Reading quality scores. */
    } state = FASTQ_STATE_ID1;

    seq->id1.n = seq->seq.n = seq->id2.n = seq->qual.n = 0;
    char* end = f->buf + f->readlen;
    do {
        while (f->next < end) {
            /* Consume pointless special characters prefixing IDs */
            if ((state == FASTQ_STATE_ID1 && f->linestart && f->next[0] == '@') ||
                (state == FASTQ_STATE_ID2 && f->linestart && f->next[0] == '+')) {
                f->linestart = false;
                ++f->next;
                continue;
            }

            char* u = memchr(f->next, '\n', end - f->next);
            if (u == NULL) {
                f->linestart = false;
                u = end;
            }
            else f->linestart = true;

            switch (state) {
                case FASTQ_STATE_ID1:
                    str_append(&seq->id1, f->next, u - f->next);
                    if (f->linestart) state = FASTQ_STATE_SEQ;
                    break;

                case FASTQ_STATE_SEQ:
                    str_append(&seq->seq, f->next, u - f->next);
                    if (f->linestart) state = FASTQ_STATE_ID2;
                    break;

                case FASTQ_STATE_ID2:
                    str_append(&seq->id2, f->next, u - f->next);
                    if (f->linestart) state = FASTQ_STATE_QUAL;
                    break;

                case FASTQ_STATE_QUAL:
                    str_append(&seq->qual, f->next, u - f->next);
                    if (f->linestart) {
                        f->next = u + 1;
                        return true;
                    }
                    break;
            }

            f->next = u + 1;
        }

        /* Try to read more. */
        f->readlen = fread(f->buf, 1, parser_buf_size, f->file);
        f->next = f->buf;
        end = f->buf + f->readlen;
    } while (f->readlen);

    return false;
}


void fastq_rewind(fastq_t* f)
{
    rewind(f->file);
    f->next = f->buf;
    f->readlen = 0;
}


void fastq_print(FILE* fout, const seq_t* seq)
{
    fprintf(fout, "@%s\n%s\n+%s\n%s\n",
                  seq->id1.s,
                  seq->seq.s,
                  seq->id2.s,
                  seq->qual.s );
}


