
#include <string.h>

#include "kmerset.h"
#include "misc.h"

/* Prime numbers that a near powers of two, suitable for hash table sizes, when
 * using quadratic probing. */
#define NUM_PRIMES 28
static const uint32_t primes[NUM_PRIMES] = {
           53U,         97U,        193U,        389U,
          769U,       1543U,       3079U,       6151U,
        12289U,      24593U,      49157U,      98317U,
       196613U,     393241U,     786433U,    1572869U,
      3145739U,    6291469U,   12582917U,   25165843U,
     50331653U,  100663319U,  201326611U,  402653189U,
    805306457U, 1610612741U, 3221225473U, 4294967291U };


/* Load factor before resize. */
static const double MAX_LOAD  = 0.7;


/* simple quadratic probing */
static uint32_t probe(uint32_t h, uint32_t i)
{
    static const uint32_t c1 = 2;
    static const uint32_t c2 = 2;

    return h + i/c1 + (i*i)/c2;
}


typedef struct kmerset_cell_t_
{
    kmer_t x;
    uint32_t idx;
} kmerset_cell_t;


struct kmerset_t_
{
    kmerset_cell_t* xs;

    /* Size of xs. */
    size_t size;

    /* Number of non-empty cells. */
    size_t n;

    /* Maxmimum load before resizing. */
    size_t max_n;
};


kmerset_t* kmerset_alloc()
{
    kmerset_t* H = malloc_or_die(sizeof(kmerset_t));
    H->n = 0;
    H->size = 0;
    H->max_n = (uint32_t) (MAX_LOAD * (double) primes[H->size]);
    H->xs = malloc_or_die(primes[H->size] * sizeof(kmerset_cell_t));
    memset(H->xs, 0, primes[H->size] * sizeof(kmerset_cell_t));
    return H;
}


void kmerset_free(kmerset_t* H)
{
    free(H->xs);
    free(H);
}


size_t kmerset_size(const kmerset_t* H)
{
    return H->n;
}


static void kmerset_expand(kmerset_t* H)
{
    ++H->size;
    kmerset_cell_t* xs = malloc_or_die(primes[H->size] * sizeof(kmerset_cell_t));
    memset(xs, 0, primes[H->size] * sizeof(kmerset_cell_t));

    uint32_t probe_num, h, k;
    size_t i;
    for (i = 0; i < primes[H->size - 1]; ++i) {
        if (H->xs[i].idx != 0) {
            probe_num = 1;
            h = kmer_hash(H->xs[i].x);
            k = h % primes[H->size];
            while (true) {
                if (xs[k].idx == 0) {
                    xs[k].x   = H->xs[i].x;
                    xs[k].idx = H->xs[i].idx;
                    break;
                }

                k = probe(h, ++probe_num) % primes[H->size];
            }
        }
    }

    free(H->xs);
    H->xs = xs;
    H->max_n = (uint32_t) (MAX_LOAD * (double) primes[H->size]);
}


void kmerset_add(kmerset_t* H, kmer_t x)
{
    if (H->n >= H->max_n) {
        kmerset_expand(H);
    }

    uint32_t probe_num = 1;
    uint32_t h = kmer_hash(x);
    uint32_t k = h % primes[H->size];

    while (true) {
        if (H->xs[k].idx == 0) {
            H->xs[k].x = x;
            H->xs[k].idx = ++H->n;
            return;
        }
        else if (H->xs[k].x == x) {
            return;
        }

        k = probe(h, ++probe_num) % primes[H->size];
    }
}


uint32_t kmerset_get(const kmerset_t* H, kmer_t x)
{
    uint32_t probe_num = 1;
    uint32_t h = kmer_hash(x);
    uint32_t k = h % primes[H->size];

    while (true) {
        if (H->xs[k].idx == 0) {
            return 0;
        }
        else if (H->xs[k].x == x) {
            return H->xs[k].idx;
        }

        k = probe(h, ++probe_num) % primes[H->size];
        if (k == h % primes[H->size]) return 0;
    }
}


