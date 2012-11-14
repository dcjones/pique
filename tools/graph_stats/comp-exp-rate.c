
/* Compute the "component expansion rate" (for lack of a better term).
 *
 * High weight edges are repeatedly removed from the graph and the number of
 * connected components is counted.
 *
 * The rate at which the number of components increase gives some hint at the
 * complexity or diversity present in the sample.
 */

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static void print_usage(FILE* f)
{
        fprintf(f,
"Usage: comp-exp-rate adjmat.mm\n"
"Estimate the component expansion rate given a adjacency matrix in mm format.\n");
}


static void* malloc_or_die(size_t n)
{
    void* p = malloc(n);
    if (p == NULL) {
        fprintf(stderr, "Can't allocate %zu bytes.\n", n);
        exit(EXIT_FAILURE);
    }
    return p;
}


static void* realloc_or_die(void* p, size_t n)
{
    if ((p = realloc(p, n)) == NULL) {
        fprintf(stderr, "Can't allocate %zu bytes.\n", n);
        exit(EXIT_FAILURE);
    }
    return p;
}

typedef struct
{
    unsigned int u, /* source */
                 v, /* dest */
                 w; /* weight */
} edge_t;


static void edgecpy(edge_t* a, const edge_t* b)
{
    a->u = b->u;
    a->v = b->v;
    a->w = b->w;
}


typedef struct
{
    edge_t* es;  /* edges */
    size_t m;    /* number of edges */
    size_t size; /* reserved space */
} edge_array_t;



static void push_edge(edge_array_t* E, const edge_t* e)
{
    if (E->m >= E->size) {
        E->size *= 2;
        E->es = realloc_or_die(E->es, E->size * sizeof(edge_t));
    }

    edgecpy(&E->es[E->m++], e);
}


static int edge_cmp(const void* e1_, const void* e2_) {
    return (int) ((edge_t*) e1_)->w -
           (int) ((edge_t*) e2_)->w;
}


/* Disjoint set data structure operations implemented on a flat array ds in
 * which ds[i] points to the index of the parent node, and ds[i] = i indicates a
 * root node. */

/* Disjoint set find with path compression. */
static unsigned int disjset_find(unsigned int* ds, unsigned int i)
{
    if (ds[i] == i) {
        return i;
    } else {
        return ds[i] = disjset_find(ds, ds[i]);
    }
}


/* Disjoint set union. */
static void disjset_union(unsigned int* ds, unsigned int i, unsigned int j)
{
    unsigned int a = disjset_find(ds, i);
    unsigned int b = disjset_find(ds, j);
    ds[b] = a;
}


static int uintcmp(const void* a_, const void* b_)
{
    unsigned int a = *(unsigned int*) a_;
    unsigned int b = *(unsigned int*) b_;
    if      (a < b) return -1;
    else if (a > b) return  1;
    else            return  0;
}


/* Count connected components, where ds is a work space of size n. */
static size_t count_components(const edge_t* es, size_t n, size_t m,
                               unsigned int* ds)
{
    size_t i;
    for (i = 0; i < n; ++i) ds[i] = i;
    for (i = 0; i < m; ++i) {
        disjset_union(ds, es[i].u, es[i].v);
    }

    for (i = 0; i < n; ++i) disjset_find(ds, i);

    qsort(ds, n, sizeof(unsigned int), uintcmp);

    /* count components */
    size_t cnt = 1;
    for (i = 1; i < n; ++i) {
        if (ds[i] != ds[i - 1]) ++cnt;
    }

    return cnt;
}



int main(int argc, char* argv[])
{
    while (1) {
        int opt = getopt(argc, argv, "h");
        if (opt == -1) break;

        switch (opt) {
            case 'h':
                print_usage(stdout);
                return EXIT_SUCCESS;

            case '?':
                return EXIT_FAILURE;

            default:
                abort();
        }
    }

    if (optind >= argc) {
        print_usage(stderr);
        return EXIT_FAILURE;
    }

    const char* fn = argv[optind];
    FILE* f = fopen(fn, "r");
    if (f == NULL) {
        fprintf(stderr, "Can't open %s for reading.\n", fn);
        return EXIT_FAILURE;
    }

    fprintf(stderr, "Reading adjacency matrix ... ");

    /* some rather brittle parsing of matrix market files */
    char buffer[512];
    fgets(buffer, sizeof(buffer), f);
    if (strcmp(buffer, "%%MatrixMarket matrix coordinate integer general\n") != 0) {
        fprintf(stderr, "Error: Incorrectly formatted matrix market file.\n");
        return EXIT_FAILURE;
    }

    unsigned int n, n2, m;
    fscanf(f, "%u %u %u\n", &n, &n2, &m);

    edge_array_t E;
    E.size = m;
    E.m = 0;
    E.es = malloc_or_die(E.size * sizeof(edge_t));
    edge_t e;

    while (fgets(buffer, sizeof(buffer), f)) {
        sscanf(buffer, "%u %u %u\n", &e.u, &e.v, &e.w);
        e.u -= 1; e.v -= 1; /* make 0-based */
        push_edge(&E, &e);
    }
    fclose(f);
    qsort(E.es, E.m, sizeof(edge_t), edge_cmp);
    fprintf(stderr, "done. (%zu edges)\n", E.m);

    unsigned int* ds = malloc_or_die(n * sizeof(unsigned int));
    size_t i;
    for (i = 0; i < E.m; ++i) {
        size_t c = count_components(E.es + i, n, E.m - i, ds);
        printf("%zu\n", c);
    }


    free(ds);
    free(E.es);

    return EXIT_SUCCESS;
}






