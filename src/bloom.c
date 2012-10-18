
#include <assert.h>
#include <string.h>
#include <pthread.h>

#include "bloom.h"
#include "misc.h"


/* the number of subtables, hard-coded so I can use stack space in a few places
 * */
#define NUM_SUBTABLES 4

/* careful, these numbers should not be changed independent of each other */
static const size_t   fingerprint_bits = 14;
static const uint32_t fingerprint_mask = 0xfffc00;
static const size_t   counter_bits     = 10;
static const uint32_t counter_mask     = 0x0003ff;
static const size_t   cell_bytes       = 3;
static const size_t   blocks_per_lock  = 16;


static uint32_t get_cell_count(uint8_t* c)
{
    return (*(uint32_t*) c) & counter_mask;
}


static void set_cell_count(uint8_t* c, uint32_t cnt)
{
    (*(uint32_t*) c) = ((*(uint32_t*) c) & fingerprint_mask) | (cnt & counter_mask);
}


struct bloom_t_
{
    /* pointers into T, to save a little computation */
    uint8_t* subtables[NUM_SUBTABLES];

    /* Mutexes each locking a group of blocks_per_lock */
    pthread_mutex_t* mutexes[NUM_SUBTABLES];

    /* number of buckets per subtable */
    size_t n;

    /* number of cells per bucket */
    size_t m;
};



bloom_t* bloom_alloc(size_t n, size_t m)
{
    bloom_t* B = malloc(sizeof(bloom_t));
    B->n = n;
    B->m = m;

    size_t subtable_size = n * m * cell_bytes;
    /* this is ceiling(n * m / blocks_per_lock) */
    size_t mutex_count = (n * m + blocks_per_lock - 1) / blocks_per_lock;

    size_t i, j;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        B->subtables[i] = malloc_or_die(subtable_size);
        memset(B->subtables[i], 0, subtable_size);

        B->mutexes[i] = malloc_or_die(mutex_count * sizeof(pthread_mutex_t));
        for (j = 0; j < mutex_count; ++j) {
            pthread_mutex_init_or_die(&B->mutexes[i][j], NULL);
        }
    }

    return B;
}


bloom_t* bloom_copy(const bloom_t* B)
{
    bloom_t* C = malloc_or_die(sizeof(bloom_t));
    C->n = B->n;
    C->m = B->m;

    size_t subtable_size = C->n * C->m * cell_bytes;
    size_t mutex_count = (C->n * C->m + blocks_per_lock - 1) / blocks_per_lock;

    size_t i, j;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        C->subtables[i] = malloc_or_die(subtable_size);
        memcpy(C->subtables[i], B->subtables[i], subtable_size);

        C->mutexes[i] = malloc_or_die(mutex_count * sizeof(pthread_mutex_t));
        for (j = 0; j < mutex_count; ++j) {
            pthread_mutex_init_or_die(&C->mutexes[i][j], NULL);
        }
    }

    return C;
}


void bloom_clear(bloom_t* B)
{
    size_t i;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        memset(B->subtables[i], 0, B->n * B->m * cell_bytes);
    }
}


void bloom_free(bloom_t* B)
{
    if (B == NULL) return;

    size_t mutex_count = (B->n * B->m + blocks_per_lock - 1) / blocks_per_lock;
    size_t i, j;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        free(B->subtables[i]);

        for (j = 0; j < mutex_count; ++j) {
            pthread_mutex_destroy(&B->mutexes[i][j]);
        }
        free(B->mutexes[i]);
    }

    free(B);
}


/* Find the subtable i and cell j containing the given key x.
 *
 * Args:
 *   B: A bloom fliter.
 *   x: The kmer to finde.
 *   i_: If located, the index of the subtable is output here.
 *   j_: If located, the offset into the subtable is output here.
 *   locket_mutex: If the find was successful, this is set to a locked mutex
 *                  protecting the location of the key. The caller is
 *                  responsible for unlocking.
 * Returns:
 *   true if the key was found.
 */
static bool bloom_find(const bloom_t* B, kmer_t x,
                       size_t* i_, size_t* j_,
                       pthread_mutex_t** locked_mutex)
{
    const size_t bytes_per_bucket = B->m * cell_bytes;

    uint64_t h1, h0 = kmer_hash(x);
    uint32_t fp = h0 & (uint64_t) fingerprint_mask;
    uint64_t hs[NUM_SUBTABLES];

    h1 = h0;
    size_t i;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        h1 = hs[i] = kmer_hash_mix(h0, h1);
        hs[i] %= B->n;
        prefetch(&B->subtables[i][hs[i] * bytes_per_bucket], 0, 0);
        prefetch(&B->mutexes[i][hs[i] / blocks_per_lock], 0, 0);
    }

    uint8_t *c, *c_start, *c_end;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        *locked_mutex = &B->mutexes[i][hs[i] / blocks_per_lock];
        pthread_mutex_lock(*locked_mutex);

        c = c_start = &B->subtables[i][hs[i] * bytes_per_bucket];
        c_end = c_start + bytes_per_bucket;

        /* scan through cells */
        while (c < c_end) {
            if (((*(uint32_t*) c) & fingerprint_mask) == fp) {
                *i_ = i;
                *j_ = c - B->subtables[i];
                return true;
            }

            c += cell_bytes;
        }

        pthread_mutex_unlock(*locked_mutex);
    }

    *locked_mutex = NULL;
    return false;
}


unsigned int bloom_get(bloom_t* B, kmer_t x)
{
    pthread_mutex_t* mutex;
    size_t i, j;
    if (bloom_find(B, x, &i, &j, &mutex)) {
        size_t count = get_cell_count(&B->subtables[i][j]);
        pthread_mutex_unlock(mutex);
        return count;
    }
    else return 0;
}


void bloom_del(bloom_t* B, kmer_t x)
{
    pthread_mutex_t* mutex;
    size_t i, j;
    if (bloom_find(B, x, &i, &j, &mutex)) {
        uint32_t* c = (uint32_t*) &B->subtables[i][j];
        *c &= ~(fingerprint_mask | counter_mask);
        pthread_mutex_unlock(mutex);
    }
}


unsigned int bloom_inc(bloom_t* B, kmer_t x)
{
    return bloom_add(B, x, 1);
}



/* Add d to the count for the key x.
 *
 * Args:
 *   B: The bloom filter.
 *   x: A key to increase.
 *   d: Delta by which to increase the key's count.
 *
 * Returns:
 *   The new count for the cell, or 0 if there was not space to place it.
 */
unsigned int bloom_add(bloom_t* B, kmer_t x, unsigned int d)
{
    /* We can't quite use bloom_find here since we have to keep track of
     * candidate cells, and more importantly, keep them locked. */

    const size_t bytes_per_bucket = B->m * cell_bytes;

    uint64_t h1, h0 = kmer_hash(x);
    uint32_t fp = h0 & (uint64_t) fingerprint_mask;
    uint64_t hs[NUM_SUBTABLES];

    h1 = h0;
    size_t i;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        h1 = hs[i] = kmer_hash_mix(h0, h1);
        hs[i] %= B->n;
        prefetch(&B->subtables[i][hs[i] * bytes_per_bucket], 0, 0);
        prefetch(&B->mutexes[i][hs[i] / blocks_per_lock], 0, 0);
    }

    uint8_t* cells[NUM_SUBTABLES];
    size_t bucket_sizes[NUM_SUBTABLES];
    pthread_mutex_t* locked_mutexes[NUM_SUBTABLES];
    memset(locked_mutexes, 0, NUM_SUBTABLES * sizeof(pthread_mutex_t*));

    uint32_t count;
    uint8_t *c, *c_start, *c_end;
    uint32_t cell_fp;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        locked_mutexes[i] = &B->mutexes[i][hs[i] / blocks_per_lock];
        pthread_mutex_lock(locked_mutexes[i]);

        c = c_start = &B->subtables[i][hs[i] * bytes_per_bucket];
        c_end = c_start + bytes_per_bucket;

        while (c < c_end) {
            cell_fp = (*(uint32_t*) c) & fingerprint_mask;

            /* Key found. */
            if (cell_fp == fp) {
                count = get_cell_count(c);
                if (count + d < counter_mask) set_cell_count(c, count + d);
                else set_cell_count(c, counter_mask);

                size_t k;
                for (k = 0; k <= i; ++k) {
                    if (locked_mutexes[i]) {
                        pthread_mutex_unlock(locked_mutexes[k]);
                    }
                }

                return count + d;
            }
            /* Candidate cell found. */
            else if (cell_fp == 0) {
                cells[i] = c;
                bucket_sizes[i] = (c - c_start) / cell_bytes;
                break;
            }

            c += cell_bytes;
        }

        /* full bucket */
        if (c == c_end) {
            cells[i] = NULL;
            bucket_sizes[i] = B->m;
            pthread_mutex_unlock(locked_mutexes[i]);
            locked_mutexes[i] = NULL;
        }
    }

    /* Find the least-full bucket, breaking ties to the left. (i.e., "d-left"
     * hashing). */
    size_t i_min = NUM_SUBTABLES;
    size_t min_bucket_size = B->m;
    for (i = 0; i < NUM_SUBTABLES && min_bucket_size > 0; ++i) {
        if (bucket_sizes[i] < min_bucket_size) {
            i_min = i;
            min_bucket_size = bucket_sizes[i];
        }
    }

    /* Insert if a suitable cell was found. */
    count = 0;
    if (i_min < NUM_SUBTABLES) {
        if (d > counter_mask) d = counter_mask;
        (*(uint32_t*) cells[i_min]) = fp | d; // figngerprint & count
        count = d;
    }

    /* Unlock the mutexes. */
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        if (locked_mutexes[i]) pthread_mutex_unlock(locked_mutexes[i]);
    }

    return count;
}


