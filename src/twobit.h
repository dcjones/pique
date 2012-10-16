/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * twobit:
 * Nucleotide sequences encoded two bits per nucleotide.
 */

#ifndef QUIP_TWOBIT
#define QUIP_TWOBIT

#include "kmer.h"
#include <stdlib.h>
#include <stdio.h>

typedef struct twobit_t_ twobit_t;

twobit_t* twobit_alloc();
twobit_t* twobit_alloc_n(size_t n);
void      twobit_free(twobit_t*);
twobit_t* twobit_dup(const twobit_t*);
void      twobit_clear(twobit_t*);
void      twobit_reserve(twobit_t* s, size_t seqlen);
void      twobit_free_reserve(twobit_t* s);

size_t twobit_len(const twobit_t*);
void   twobit_copy(twobit_t*, const twobit_t*);
void   twobit_copy_str(twobit_t*, const char*);
void   twobit_copy_str_n(twobit_t*, const char*, size_t);
void   twobit_append(twobit_t*, const char*);
void   twobit_append_char(twobit_t*, char);
void   twobit_append_n(twobit_t*, const char*, size_t);
void   twobit_append_kmer(twobit_t*, kmer_t x, size_t k);
void   twobit_append_twobit(twobit_t*, const twobit_t*);
void   twobit_reverse(twobit_t*);


void   twobit_setc(twobit_t*, size_t i, char);
void   twobit_set(twobit_t*, size_t i, kmer_t);
kmer_t twobit_get(const twobit_t*, size_t i);
kmer_t twobit_get_kmer(const twobit_t*, size_t i, size_t k);
kmer_t twobit_get_kmer_rev(const twobit_t* s, size_t i, size_t k);
void   twobit_print(const twobit_t*, FILE*);
void   twobit_print_stdout(const twobit_t*);
int    twobit_cmp(const twobit_t*, const twobit_t*);
void   twobit_revcomp(twobit_t* dest, const twobit_t* src);

uint32_t twobit_hash(const twobit_t*);
uint64_t twobit_crc64_update(const twobit_t*, uint64_t crc);

/* Count mismatches (i.e. hamming distange) between a query and subject,
 * with the query places at the given offset in the subject. */
uint32_t twobit_mismatch_count(const twobit_t* subject,
                               const twobit_t* query,
                               size_t spos, uint32_t max_miss);

#endif

