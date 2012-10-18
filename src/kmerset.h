
#ifndef PIQUE_KMERHASH
#define PIQUE_KMERHASH

#include "kmer.h"

typedef struct kmerset_t_ kmerset_t;

kmerset_t* kmerset_alloc();
void kmerset_free(kmerset_t*);

size_t kmerset_size(const kmerset_t* H);

void kmerset_add(kmerset_t* H, kmer_t x);

/* Return the (one-based) index of the kmer in the set.
 *
 * Zero is returned if the kmer is not present in the set. */
uint32_t kmerset_get(const kmerset_t* H, kmer_t x);

#endif

