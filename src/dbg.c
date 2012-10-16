
#include "bloom.h"
#include "dbg.h"
#include "kmercache.h"
#include "misc.h"

/* I'm fixing cells per block. It's not obvious the effect of changing it, so I
 * don't want to expose it as an option. */
static const size_t cells_per_bucket = 8;

/* Maximum number of seeds we might accumulated. */
static const size_t max_seeds = 250000;

struct dbg_t_
{
    /* Bloom filter to accumalate k-mer statistics. */
    bloom_t* B;

    /* k-mer size */
    size_t k;

    /* k-mer mask */
    kmer_t mask;

    /* A leaky hash table of k-mer seeds used as starting points for traversing
     * the graph. */
    kmercache_t* seeds;
};


dbg_t* dbg_alloc(size_t n, size_t k)
{
    dbg_t* G = malloc_or_die(sizeof(dbg_t));

    size_t num_buckets = n / 4 / cells_per_bucket; /* assuming 4 subtables. */
    G->B = bloom_alloc(num_buckets , cells_per_bucket);
    G->k = k;
    G->mask = (1 << (2 * k)) - 1;
    G->seeds = kmercache_alloc(max_seeds);
    return G;
}


void dbg_free(dbg_t* G)
{
    bloom_free(G->B);
    kmercache_free(G->seeds);
    free(G);
}


void dbg_add_twobit_seq(dbg_t* G, const twobit_t* seq)
{
    size_t i, len = twobit_len(seq);
    kmer_t x = 0, y;
    for (i = 0; i < len; ++i) {
        x = ((x << 2) | twobit_get(seq, i)) & G->mask;

        if (i + 1 >= G->k) {
            y = kmer_canonical(x, G->k);
            bloom_add(G->B, y, 1);
            kmercache_inc(G->seeds, y);
        }
    }
}


void dbg_dump(const dbg_t* G, FILE* fout)
{
    /* What are we actually doing here?
     *
     * We have to use seeded traversals sirce the dlcbf doesn't allow us to just
     * list k-mers that are present.
     *
     * So be it, let's figure out a suitable seeding strategy.
     *
     * In quip, we just use the last n reads, where n is 200,000.
     *
     * That's a little shoddy, let's try to come up with something more
     * interesting.
     *
     * Here's what I'm thinking: we use a fixed size hash table that records
     * counts. Collisions are resolved probabalisticlly: if the count is high,
     * it has a increasingly small probability of getting booted out.
     * */
}


