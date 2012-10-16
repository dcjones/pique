/*
 * This file is part of pique.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/*
 * debruijn:
 * A probabalistic de bruijn graph implemented with a dlcbf.
 */

#ifndef PIQUE_DEBRUIJN
#define PIQUE_DEBRUIJN

#include "twobit.h"

typedef struct dbg_t_ dbg_t;

/* A allocate a De Bruijn graph.
 *
 * Args:
 *   n: Reserve space for this many unique k-mers.
 *
 * Returns:
 *   An allocated graph.
 */
dbg_t* dbg_alloc(size_t n, size_t k);


/* Free a graph allocated with dbg_alloc. */
void dbg_free(dbg_t* G);


/* */
void dbg_add_twobit_seq(dbg_t* G, const twobit_t* seq);


#endif

