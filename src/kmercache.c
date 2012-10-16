
#include <math.h>
#include <string.h>

#include "kmercache.h"
#include "misc.h"

/* Coarseness of the locking. */
static const size_t cells_per_mutex = 16;

/* Base probability of the current occupant being booted upon collision. */
static const double base_rep_pr = 0.9;


kmercache_t* kmercache_alloc(size_t n)
{
    kmercache_t* C = malloc_or_die(sizeof(kmercache_t));
    C->n = n;
    C->xs = malloc_or_die(n * sizeof(kmercache_cell_t));
    memset(C->xs, 0, n * sizeof(kmercache_cell_t));

    size_t num_mutexes = (n + cells_per_mutex - 1) / cells_per_mutex;
    C->mutexes = malloc_or_die(num_mutexes * sizeof(pthread_mutex_t));
    size_t i;
    for (i = 0; i < num_mutexes; ++i) {
        pthread_mutex_init_or_die(&C->mutexes[i], NULL);
    }

    C->rng = rng_alloc(12345);
    pthread_mutex_init_or_die(&C->rng_mutex, NULL);

    return C;
}


void kmercache_free(kmercache_t* C)
{
    size_t i, num_mutexes = (C->n + cells_per_mutex - 1) / cells_per_mutex;
    for (i = 0; i < num_mutexes; ++i) {
        pthread_mutex_destroy(&C->mutexes[i]);
    }
    free(C->mutexes);
    free(C->xs);
    rng_free(C->rng);
    pthread_mutex_destroy(&C->rng_mutex);
    free(C);
}


uint32_t kmercache_inc(kmercache_t* C, kmer_t x)
{
    uint64_t i = kmer_hash(x) % C->n;
    uint32_t count = 0;
    pthread_mutex_lock(&C->mutexes[i / cells_per_mutex]);

    if (C->xs[i].x == x) {
        if (C->xs[i].count < UINT32_MAX) ++C->xs[i].count;
        count = C->xs[i].count;
    }
    else {
        double pr = pow(base_rep_pr, (double) C->xs[i].count);
        pthread_mutex_lock(&C->rng_mutex);
        double r = rng_get_double(C->rng);
        pthread_mutex_unlock(&C->rng_mutex);
        if (r < pr) {
            C->xs[i].x = x;
            count = C->xs[i].count = 1;
        }
    }

    pthread_mutex_unlock(&C->mutexes[i / cells_per_mutex]);
    return count;
}

