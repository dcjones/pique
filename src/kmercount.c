
#include <string.h>

#include "kmercount.h"
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


/* simple quadratic probing */
static uint32_t probe(uint32_t h, uint32_t i)
{
    static const uint32_t c1 = 2;
    static const uint32_t c2 = 2;

    return h + i/c1 + (i*i)/c2;
}


/* Load factor before resize. */
static const double MAX_LOAD  = 0.7;


static kmercount_table_t* kmercount_table_alloc()
{
    kmercount_table_t* T = malloc_or_die(sizeof(kmercount_table_t));
    T->n = 0;
    T->size = 0;
    T->max_n = (uint32_t) (MAX_LOAD * (double) primes[T->size]);
    T->xs = malloc_or_die(primes[T->size] * sizeof(kmercount_cell_t));
    memset(T->xs, 0, primes[T->size] * sizeof(kmercount_cell_t));
    pthread_mutex_init_or_die(&T->mut, NULL);
    return T;
}


static void kmercount_table_free(kmercount_table_t* T)
{
    if (T) {
        free(T->xs);
        pthread_mutex_destroy(&T->mut);
        free(T);
    }
}


static void kmercount_table_expand(kmercount_table_t* T)
{
    ++T->size;
    kmercount_cell_t* xs = malloc_or_die(primes[T->size] * sizeof(kmercount_cell_t));
    memset(xs, 0, primes[T->size] * sizeof(kmercount_cell_t));

    uint32_t probe_num, h, k;
    size_t i;
    for (i = 0; i < primes[T->size - 1]; ++i) {
        if (T->xs[i].count != 0) {
            probe_num = 1;
            h = kmer_hash(T->xs[i].x);
            k = h % primes[T->size];
            while (true) {
                if (xs[k].count == 0) {
                    xs[k].x      = T->xs[i].x;
                    xs[k].count = T->xs[i].count;
                    break;
                }

                k = probe(h, ++probe_num) % primes[T->size];
            }
        }
    }

    free(T->xs);
    T->xs = xs;
    T->max_n = (uint32_t) (MAX_LOAD * (double) primes[T->size]);
}


static void kmercount_table_add(kmercount_table_t* T, kmer_t x,
                                uint32_t delta, uint32_t h)
{
    pthread_mutex_lock(&T->mut);

    if (T->n >= T->max_n) {
        kmercount_table_expand(T);
    }

    uint32_t probe_num = 1;
    uint32_t k = h % primes[T->size];
    while (true) {
        if (T->xs[k].count == 0) {
            T->xs[k].x = x;
            T->xs[k].count = delta;
            ++T->n;
            break;
        }
        else if (T->xs[k].x == x) {
            T->xs[k].count += delta;
            break;
        }

        k = probe(h, ++probe_num) % primes[T->size];
    }

    pthread_mutex_unlock(&T->mut);
}


static void kmercount_table_set(kmercount_table_t* T, kmer_t x,
                                uint32_t count, uint32_t h)
{
    pthread_mutex_lock(&T->mut);

    if (T->n >= T->max_n) {
        kmercount_table_expand(T);
    }

    uint32_t probe_num = 1;
    uint32_t k = h % primes[T->size];
    while (true) {
        if (T->xs[k].count == 0) {
            T->xs[k].x = x;
            T->xs[k].count = count;
            ++T->n;
            break;
        }
        else if (T->xs[k].x == x) {
            T->xs[k].count = count;
            break;
        }

        k = probe(h, ++probe_num) % primes[T->size];
    }

    pthread_mutex_unlock(&T->mut);
}


static uint32_t kmercount_table_get(kmercount_table_t* T, kmer_t x, uint32_t h)
{
    pthread_mutex_lock(&T->mut);

    uint32_t probe_num = 1;
    uint32_t k = h % primes[T->size];
    uint32_t count = 0;

    while (true) {
        if (T->xs[k].count == 0) {
            break;
        }
        else if (T->xs[k].x == x) {
            count = T->xs[k].count;
            break;
        }

        k = probe(h, ++probe_num) % primes[T->size];
    }

    pthread_mutex_unlock(&T->mut);

    return count;
}


kmercount_t* kmercount_alloc(size_t k)
{
    kmercount_t* C = malloc_or_die(sizeof(kmercount_t));
    C->k = k;
    C->subtables = malloc_or_die(k * sizeof(kmercount_table_t*));
    size_t i;
    for (i = 0; i < k; ++i) {
        C->subtables[i] = kmercount_table_alloc();
    }
    return C;
}


void kmercount_free(kmercount_t* C)
{
    if (C) {
        size_t i;
        for (i = 0; i < C->k; ++i) {
            kmercount_table_free(C->subtables[i]);
        }
        free(C->subtables);
        free(C);
    }
}


void kmercount_add(kmercount_t* C, kmer_t x, uint32_t delta)
{
    uint64_t h = kmer_hash(x);

    /* totally ad-hoc allocation of k-mers to subtables */
    uint64_t i = (h >> 32) % C->k;

    kmercount_table_add(C->subtables[i], x, delta, (uint32_t) h);
}


uint32_t kmercount_get(kmercount_t* C, kmer_t x)
{
    uint64_t h = kmer_hash(x);

    /* totally ad-hoc allocation of k-mers to subtables */
    uint64_t i = (h >> 32) % C->k;

    return kmercount_table_get(C->subtables[i], x, (uint32_t) h);
}


void kmercount_set(kmercount_t* C, kmer_t x, uint32_t count)
{
    uint64_t h = kmer_hash(x);

    /* totally ad-hoc allocation of k-mers to subtables */
    uint64_t i = (h >> 32) % C->k;

    kmercount_table_set(C->subtables[i], x, count, (uint32_t) h);
}


