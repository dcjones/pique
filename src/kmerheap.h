
/*
 * This file is part of pique
 *
 * Copyright (c) 2013 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/* kmerheap:
 * A simple fixed-size heap used to find the top m kmers.
 * */

#ifndef PIQUE_KMERHEAP_H
#define PIQUE_KMERHEAP_H

#include "kmercount.h"

typedef struct kmerheap_t_
{
    kmercount_cell_t* xs;

    /* size of the heap */
    size_t n;

    /* number of non-zero cells */
    size_t m;

} kmerheap_t;


kmerheap_t* kmerheap_alloc(size_t n);
void kmerheap_free(kmerheap_t*);
void kmerheap_add(kmerheap_t* H, kmer_t x, uint32_t count);
bool kmerheap_pop(kmerheap_t* H, kmer_t* x, uint32_t* count);

#endif

