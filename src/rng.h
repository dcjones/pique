/*
 * This file is part of pique.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/* rng: A simple pseudo random number generator. */

#ifndef PIQUE_RNG_H
#define PIQUE_RNG_H

#include <stdint.h>

typedef struct rng_t_ rng_t;

rng_t* rng_alloc(uint32_t seed);
void rng_free(rng_t* rng);

/* Get a random uint32_t in [0, UINT32_MAX]. */
uint32_t rng_get(rng_t* rng);

/* Get a random double in [0, 1]. */
double rng_get_double(rng_t* rng);

#endif

