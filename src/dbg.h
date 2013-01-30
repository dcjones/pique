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

#ifndef PIQUE_DBG
#define PIQUE_DBG

#include "twobit.h"
#include "rng.h"

typedef struct dbg_t_ dbg_t;

/* A allocate a De Bruijn graph.
 *
 * Args:
 *   k: k-mer size
 *
 * Returns:
 *   An allocated graph.
 */
dbg_t* dbg_alloc(size_t k);


/* Free a graph allocated with dbg_alloc. */
void dbg_free(dbg_t* G);


/* Add the k-mers contained in a sequence to the de bruijn graph. */
void dbg_add_twobit_seq(dbg_t* G, const twobit_t* seq);


/* Dump the graph to a readable file. */
typedef enum {
    ADJ_GRAPH_FMT_MM,
    ADJ_GRAPH_FMT_HB
} adj_graph_fmt_t;

void dbg_dump(const dbg_t* G, FILE* fout, size_t num_threads,
              adj_graph_fmt_t fmt);


#endif

