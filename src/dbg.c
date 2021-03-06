
#include <assert.h>
#include <inttypes.h>
#include <string.h>

#include "bloom.h"
#include "dbg.h"
#include "kmercache.h"
#include "kmerset.h"
#include "misc.h"


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
    kmer_t u;
    kmer_t v;
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

    S->es[S->n].u = e->u;
    S->es[S->n].v = e->v;
    S->es[S->n].count = e->count;
    ++S->n;
}


#if 0
static bool edgestack_pop(edgestack_t* S, edge_t* e)
{
    if (S->n == 0) return false;

    --S->n;
    e->u = S->es[S->n].u;
    e->v = S->es[S->n].v;
    e->count = S->es[S->n].count;

    return true;
}
#endif


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
    G->mask = kmer_mask(k);
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


/* One thread traversing the graph. */
typedef struct dbg_dump_thread_ctx_t_
{
    bloom_t* B;
    kmerstack_t* seeds;
    size_t k;
} dbg_dump_thread_ctx_t;


/* A helper function used by dbg_dump_thread.
 *
 * Find all out-edges from the given k-mer, push to edges, and push discovered
 * nodes to S.
 * */
static void enumerate_out_edges(kmer_t u, size_t k,
                                bloom_t* B,
                                kmerstack_t* S,
                                edgestack_t* edges)
{
    kmer_t mask = kmer_mask(k);
    uint32_t count;
    edge_t e;
    e.u = u;
    kmer_t v, vc, x;
    for (x = 0; x < 4; ++x) {
        v = ((u << 2) | x) & mask;
        vc = kmer_canonical(v, k);
        count = bloom_get(B, vc);
        if (count > 0) {
            e.v = v;
            e.count = count;
            edgestack_push(edges, &e);
            kmerstack_push(S, vc);
        }
    }
}


/* A helper function used by dbg_dump_thread.
 *
 * Find all in-edges from the given k-mer, push to edges, and push discovered
 * nodes to S.
 * */
static void enumerate_in_edges(kmer_t v, size_t k, uint32_t v_count,
                               bloom_t* B, kmerstack_t* S, edgestack_t* edges)
{
    kmer_t mask = kmer_mask(k);
    uint32_t u_count;
    edge_t e;
    e.v = v;
    e.count = v_count;
    kmer_t u, uc, x;
    for (x = 0; x < 4; ++x) {;
        u = ((v >> 2) | (x << (2*(k-1)))) & mask;
        uc = kmer_canonical(u, k);
        u_count = bloom_get(B, uc);
        if (u_count > 0) {
            e.u = u;
            edgestack_push(edges, &e);
            kmerstack_push(S, uc);
        }
    }
}


/* A de bruijn graph traversal thread.
 *
 * Eeach thread starts from a seed and performs (essentially) depth-first
 * traversal, deleting nodes as it goes and pushing edges onto a stack. */
static void* dbg_dump_thread(void* arg)
{
    dbg_dump_thread_ctx_t* ctx = (dbg_dump_thread_ctx_t*) arg;
    edgestack_t* edges = edgestack_alloc();
    kmerstack_t* S = kmerstack_alloc();

    uint32_t u_count;
    kmer_t u, u_rc;
    while (kmerstack_pop(ctx->seeds, &u)) {
        u = kmer_canonical(u, ctx->k);
        do {
            u_count = bloom_get(ctx->B, u);
            if (u_count == 0) continue;

            u_rc = kmer_revcomp(u, ctx->k);

            /* TODO: It's possible here to push the same edge twice.
             * Is this ever a problem? */

            enumerate_out_edges(u, ctx->k, ctx->B, S, edges);
            enumerate_out_edges(u_rc, ctx->k, ctx->B, S, edges);

            enumerate_in_edges(u, ctx->k, u_count, ctx->B, S, edges);
            enumerate_in_edges(u_rc, ctx->k, u_count, ctx->B, S, edges);

            bloom_del(ctx->B, u);
        } while (kmerstack_pop(S, &u));
    }

    kmerstack_free(S);
    return edges;
}


/* Write a sparse adjacency matrix in matrix market exchange format. */
static void write_sparse_mm(FILE* fout,
                            size_t node_count,
                            size_t edge_count,
                            const kmerset_t* H,
                            edgestack_t* const* edges,
                            size_t num_threads)
{
    fputs("%%MatrixMarket matrix coordinate integer general\n", fout);
    fprintf(fout, "%zu %zu %zu\n", node_count, node_count, edge_count);
    unsigned int i, j, u_idx, v_idx;
    for (i = 0; i < num_threads; ++i) {
        for (j = 0; j < edges[i]->n; ++j) {
            u_idx = kmerset_get(H, edges[i]->es[j].u);
            v_idx = kmerset_get(H, edges[i]->es[j].v);
            assert(u_idx > 0);
            assert(v_idx > 0);
            fprintf(fout, "%u %u %"PRIu16"\n",
                    u_idx, v_idx, edges[i]->es[j].count);
        }
    }
}


typedef struct edgepair_t_
{
    uint32_t u, v;
    int count;
} edgepair_t;


static int edgepair_cmp(const void* a_, const void* b_)
{
    const edgepair_t* a = (edgepair_t*) a_;
    const edgepair_t* b = (edgepair_t*) b_;

    if (a->v != b->v) {
        return (long) a->v - (long) b->v;
    }
    else {
        return (long) a->u - (long) b->u;
    }
}


/* Write a sparse adjacency matrix in harwell-boeing format. */
static void write_sparse_hb(FILE* fout,
                            size_t node_count,
                            size_t edge_count,
                            const kmerset_t* H,
                            edgestack_t* const* edges,
                            size_t num_threads)
{
    /* We need to make a big array of (i,j) edge pairs then sort by j to output
     * in column format. */
    edgepair_t* pairs = malloc_or_die(edge_count * sizeof(edgepair_t));
    size_t k = 0;
    size_t i, j;
    for (i = 0; i < num_threads; ++i) {
        for (j = 0; j < edges[i]->n; ++j, ++k) {
            pairs[k].u = kmerset_get(H, edges[i]->es[j].u);
            pairs[k].v = kmerset_get(H, edges[i]->es[j].v);
            pairs[k].count = edges[i]->es[j].count;
        }
    }

    qsort(pairs, edge_count, sizeof(edgepair_t), edgepair_cmp);

    fputs("pique generated de bruijn graph adjacency matrix                        padjmat \n", fout);
    fprintf(fout, "%14zu%14zu%14zu%14zu%14zu\n",
            node_count + edge_count + edge_count + 1,
            node_count + 1, edge_count, edge_count, (size_t) 0);
    fprintf(fout, "RUA%25zu%14zu%14zu%14zu\n",
            node_count, node_count, edge_count, (size_t) 0);
    fprintf(fout, "%16s%16s%20s%20s\n", "(1I11)", "(1I11)", "(1E9.0)", "");

    /* Output pointers to columns */
    size_t col;
    for (i = 0, col = 0; i < edge_count; ++i) {
        while (pairs[i].v >= col) {
            fprintf(fout, "%11zu\n", i + 1);
            ++col;
        }
    }
    for (; col <= node_count; ++col) {
        fprintf(fout, "%11zu\n", node_count + 1);
    }

    /* Output row indexes */
    for (i = 0; i < edge_count; ++i) {
        fprintf(fout, "%11zu\n", (size_t) pairs[i].u);
    }


    /* Output data */
    for (i = 0; i < edge_count; ++i) {
        fprintf(fout, "%9d\n", pairs[i].count);
    }

    free(pairs);
}


void dbg_dump(const dbg_t* G, FILE* fout, size_t num_threads,
              adj_graph_fmt_t fmt)
{
    /* Dump seeds and sort for best-first traversal. */
    kmercache_cell_t* seeds = malloc_or_die(G->seeds->n * sizeof(kmercache_cell_t));
    memcpy(seeds, G->seeds->xs, G->seeds->n * sizeof(kmercache_cell_t));
    qsort(seeds, G->seeds->n, sizeof(kmercache_cell_t), kmer_cache_cell_cmp);

    /* We make two traversals of the graph. On the first we count the number of
     * edges and nodes, and on the second we output edges and nodes. */

    kmerstack_t* S = kmerstack_alloc();

    size_t i;
    for (i = 0; i < G->seeds->n; ++i) {
        if (seeds[i].count > 0) kmerstack_push(S, seeds[i].x);
    }

    pthread_t* threads = malloc_or_die(num_threads * sizeof(pthread_t));
    edgestack_t** edges = malloc_or_die(num_threads * sizeof(edgestack_t*));
    dbg_dump_thread_ctx_t ctx;
    ctx.B = G->B;
    ctx.seeds = S;
    ctx.k = G->k;

    for (i = 0; i < num_threads; ++i) {
        pthread_create(&threads[i], NULL, dbg_dump_thread, (void*) &ctx);
    }

    for (i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], (void**) &edges[i]);
    }

    /* Hash k-mers present in the edge list to assign matrix indexes */
    size_t j;
    size_t edge_count = 0;
    kmerset_t* H = kmerset_alloc();
    for (i = 0; i < num_threads; ++i) {
        for (j = 0; j < edges[i]->n; ++j) {
            kmerset_add(H, edges[i]->es[j].u);
            kmerset_add(H, edges[i]->es[j].v);
        }
        edge_count += edges[i]->n;
    }

    size_t node_count = kmerset_size(H);

    if (fmt == ADJ_GRAPH_FMT_HB) {
        write_sparse_hb(fout, node_count, edge_count, H, edges, num_threads);
    }
    else if (fmt == ADJ_GRAPH_FMT_MM) {
        write_sparse_mm(fout, node_count, edge_count, H, edges, num_threads);
    }

    kmerset_free(H);
    free(edges);
    free(threads);
    kmerstack_free(S);
    free(seeds);
}


