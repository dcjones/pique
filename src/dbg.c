
#include <assert.h>
#include <string.h>

#include "bloom.h"
#include "dbg.h"
#include "kmercache.h"
#include "misc.h"
#include "mmio.h"


/* Kmer stack, used for traversals of the graph. */
typedef struct kmerstack_t_
{
    kmer_t* xs;  /* kmers */
    size_t n;    /* number of elements stored. */
    size_t size; /* size allocated */
    pthread_mutex_t mutex;
} kmerstack_t;


static kmerstack_t* kmerstack_alloc()
{
    kmerstack_t* S = malloc_or_die(sizeof(kmerstack_t));
    S->n = 0;
    S->size = 1024;
    S->xs = malloc_or_die(S->size * sizeof(kmer_t));
    pthread_mutex_init_or_die(&S->mutex, NULL);
    return S;
}


static void kmerstack_free(kmerstack_t* S)
{
    pthread_mutex_destroy(&S->mutex);
    free(S->xs);
    free(S);
}


static void kmerstack_push(kmerstack_t* S, kmer_t x)
{
    pthread_mutex_lock(&S->mutex);
    if (S->n == S->size) {
        S->size *= 2;
        S->xs = realloc_or_die(S->xs, S->size * sizeof(kmer_t));
    }

    S->xs[S->n++] = x;
    pthread_mutex_unlock(&S->mutex);
}


static bool kmerstack_pop(kmerstack_t* S, kmer_t* x)
{
    pthread_mutex_lock(&S->mutex);
    if (S->n == 0) {
        pthread_mutex_unlock(&S->mutex);
        return false;
    }

    *x = S->xs[--S->n];
    pthread_mutex_unlock(&S->mutex);
    return true;
}


/* An edge used for dbg_dump. */
typedef struct edge_t_
{
    /* i and j are indexes into dlcbf, rather than kmers, to save space.
     * In most analysis we may want to do with the graph, it doesn't really
     * matter what the kmers are. */
    uint32_t i;
    uint32_t j;
    uint16_t count;
} edge_t;



/* A stack of edges used when dumping the de bruijn graph to a sparse matrix.
 * */
typedef struct edgestack_t_
{
    edge_t* es;  /* kmers */
    size_t  n;    /* number of elements stored. */
    size_t  size; /* size allocated */
} edgestack_t;


static edgestack_t* edgestack_alloc()
{
    edgestack_t* S = malloc_or_die(sizeof(edgestack_t));
    S->n = 0;
    S->size = 1024;
    S->es = malloc_or_die(S->size * sizeof(edge_t));
    return S;
}


static void edgestack_push(edgestack_t* S, const edge_t* e)
{
    if (S->n == S->size) {
        S->size *= 2;
        S->es = realloc_or_die(S->es, S->size * sizeof(edge_t));
    }

    S->es[S->n].i = e->i;
    S->es[S->n].j = e->j;
    S->es[S->n].count = e->count;
    ++S->n;
}


static bool edgestack_pop(edgestack_t* S, edge_t* e)
{
    if (S->n == 0) return false;

    --S->n;
    e->i = S->es[S->n].i;
    e->j = S->es[S->n].j;
    e->count = S->es[S->n].count;

    return true;

}


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


void dbg_add_twobit_seq(dbg_t* G, rng_t* rng, const twobit_t* seq)
{
    size_t i, len = twobit_len(seq);
    kmer_t x = 0, y;
    for (i = 0; i < len; ++i) {
        x = ((x << 2) | twobit_get(seq, i)) & G->mask;

        if (i + 1 >= G->k) {
            y = kmer_canonical(x, G->k);
            bloom_add(G->B, y, 1);
            kmercache_inc(G->seeds, rng, y);
        }
    }
}


static int kmer_cache_cell_cmp(const void* a, const void* b)
{
    uint32_t ca = ((kmercache_cell_t*) a)->count;
    uint32_t cb = ((kmercache_cell_t*) b)->count;

    if      (ca  < cb) return -1;
    else if (ca == cb) return  0;
    else               return  1;
}


/* No, theres a big problem with this: because we are launching many threads,
 * nondeterminism becomes an issue when doing two passes. We might count n
 * edges once and n' != n another time depending on how the threads are
 * executed.
 *
 * Solution?
 * One solution is to not make two passes, but build up a list of edges.
 * Ok, I guess that's what we have to do.
 */


/* One thread traversing the graph. */
typedef struct dbg_dump_thread_ctx_t_
{
    kmerstack_t* seeds;

} dbg_dump_thread_ctx_t;

static void* dbg_dump_thread(void* arg)
{
    dbg_dump_thread_ctx_t* ctx = (dbg_dump_thread_ctx_t*) arg;
    edgestack_t* edges = edgestack_alloc();
    kmerstack_t* S = kmerstack_alloc();

    kmer_t x;
    while (kmerstack_pop(ctx->seeds, &x)) {
        /* TODO: push kmers to S for forward traversal. */
        /* TODO: push kmers to S for backwards traversal. */
    }

    kmerstack_free(S);
    return edges;
}


void dbg_dump(const dbg_t* G, FILE* fout)
{
    /* Dump seeds and sort for best-first traversal. */
    kmercache_cell_t* seeds = malloc_or_die(G->seeds->n * sizeof(kmercache_cell_t));
    memcpy(seeds, G->seeds->xs, G->seeds->n * sizeof(kmercache_cell_t));
    qsort(seeds, G->seeds->n, sizeof(kmercache_cell_t), kmer_cache_cell_cmp);

    /* We make two traversals of the graph. On the first we count the number of
     * edges and nodes, and on the second we output edges and nodes. */

    /* Pass 1 */
    size_t edge_count = 0;
    size_t node_count = 0;

    kmerstack_t* S = kmerstack_alloc();

    size_t i;
    for (i = 0; i < G->seeds->n; ++i) {
        if (seeds[i].count > 0) kmerstack_push(S, seeds[i].x);
    }

    /* Ok, now we need to start t threads, each thread is going to pop nodes
     * from the stack, query the bloom filter, and push nodes back onto the
     * stack.
     *
     * If there proves to be too much contention, each thread could maintain its
     * own stack and only draw from the global one when its stack is empty.
     *
     * Yeah, that's good.
     * */


    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_integer(&matcode);
    mm_set_symmetric(&matcode);
    mm_write_banner(fout, matcode);
    mm_write_mtx_crd_size(fout, node_count, node_count, edge_count);



    /* Pass 2 */

    for (i = 0; i < G->seeds->n; ++i) {
        if (seeds[i].count > 0) kmerstack_push(S, seeds[i].x);
    }

    /* TODO: launch traversal threads that will write edges to stdout.
     *
     */


    kmerstack_free(S);
    free(seeds);
}


