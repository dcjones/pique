/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/*
 * kmercache:
 * A probabilistic leaky hash table. The table is of a fixed size, when a
 * collision occurs, the current occupied is kicked out randomly with
 * probability that decreases when higher counts.
 */

#ifndef PIQUE_KMERCACHE_H
#define PIQUE_KMERCACHE_H

#include <pthread.h>

#include "kmer.h"
#include "rng.h"

typedef struct kmercache_cell_t_
{
    kmer_t x;
    uint32_t count;
} kmercache_cell_t;


typedef struct kmercache_t_
{
    kmercache_cell_t* xs;
    pthread_mutex_t* mutexes;
    size_t n;

    rng_t* rng;
    pthread_mutex_t rng_mutex;
} kmercache_t;


/* Allocated a kmer cache with n cells. */
kmercache_t* kmercache_alloc(size_t n);


/* Free a kmer cache allocated with kmercache_alloc. */
void kmercache_free(kmercache_t* C);


/* Increment the count of the the key x.
 *
 * The key is added if it's not present.
 *
 * Args:
 *   C: kmer cache
 *   x: key
 *
 * Returns:
 *   The new count associated with the key, which can be 0 if the key could not
 *   be inserted.
 */
uint32_t kmercache_inc(kmercache_t* C, rng_t* rng, kmer_t x);

#endif

