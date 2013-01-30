
/*
 * This file is part of pique
 *
 * Copyright (c) 2013 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/* kmercount:
 * A simple hash table counting k-mer occurances.
 * */

#ifndef PIQUE_KMERCOUNT_H
#define PIQUE_KMERCOUNT_H

#include <pthread.h>
#include "kmer.h"

typedef struct kmercount_cell_t
{
    kmer_t x;
    uint32_t count;
} kmercount_cell_t;


/* A simple hash table for k-mers. */
typedef struct kmercount_table_t_
{
    /* cells */
    kmercount_cell_t* xs;

    /* occupied cells */
    size_t n;

    /* reserved space */
    size_t size;

    /* maxmimum load before resizing. */
    size_t max_n;

    /* global lock for any access */
    pthread_mutex_t mut;
} kmercount_table_t;


/* A concurrent hash table distributed across k subtables. */
typedef struct kmercount_t_
{
    kmercount_table_t** subtables;

    /* number of subtables */
    size_t k;
} kmercount_t;


kmercount_t* kmercount_alloc(size_t k);
void kmercount_free(kmercount_t*);

/* Add delta to the current count of the kmer x, inserting it if it is not
 * already in the table. */
void kmercount_add(kmercount_t* C, kmer_t x, uint32_t delta);

void kmercount_set(kmercount_t* C, kmer_t x, uint32_t count);

uint32_t kmercount_get(kmercount_t* C, kmer_t x);


#endif

