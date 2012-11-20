
#ifndef FASTQ_TOOLS_RNG_H
#define FASTQ_TOOLS_RNG_H

typedef struct rng_t_ rng_t;

rng_t* rng_alloc();
void   rng_free(rng_t*);
void rng_seed(rng_t*, unsigned long seed);

/* Uniform integer in [0, k-1] */
unsigned long rng_uniform_int(rng_t*, unsigned long k);

#endif

