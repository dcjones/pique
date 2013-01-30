
#include <string.h>

#include "kmerheap.h"
#include "misc.h"

kmerheap_t* kmerheap_alloc(size_t n)
{
    kmerheap_t* H = malloc_or_die(sizeof(kmerheap_t));
    H->xs = malloc_or_die(n * sizeof(kmercount_cell_t));
    H->n = n;
    H->m = 0;
    memset(H->xs, 0, n * sizeof(kmercount_cell_t));
    return H;
}


void kmerheap_free(kmerheap_t* H)
{
    if (H) {
        free(H->xs);
        free(H);
    }
}


static size_t hparent(size_t i)
{
    return (i - 1) / 2;
}


static size_t hleft(size_t i)
{
    return 2 * i + 1;
}


static size_t hright(size_t i)
{
    return 2 * i + 2;
}


static void hswap(kmerheap_t* H, size_t i, size_t j)
{
    kmercount_cell_t tmp;
    tmp.x     = H->xs[i].x;
    tmp.count = H->xs[i].count;

    H->xs[i].x     = H->xs[j].x;
    H->xs[i].count = H->xs[j].count;

    H->xs[j].x     = tmp.x;
    H->xs[j].count = tmp.count;
}


void kmerheap_add(kmerheap_t* H, kmer_t x, uint32_t count)
{
    size_t i = H->m;
    H->xs[i].x     = x;
    H->xs[i].count = count;

    /* percolate up */
    size_t j;
    while (i > 0) {
        j = hparent(i);
        if (H->xs[i].count > H->xs[j].count) {
            hswap(H, i, j);
            i = j;
        }
        else break;
    }

    if (H->m < H->n - 1) ++H->m;
    else {
        /* boot out last place if we are full */
        H->xs[H->m].count = 0;
    }
}


bool kmerheap_pop(kmerheap_t* H, kmer_t* x, uint32_t* count)
{
    if (H->m == 0) return false;

    *x     = H->xs[0].x;
    *count = H->xs[0].count;

    --H->m;
    H->xs[0].x     = H->xs[H->m].x;
    H->xs[0].count = H->xs[H->m].count;

    /* percolate down */
    size_t l, r, j, i = 0;
    while (true) {
        l = hleft(i);
        r = hright(i);

        if (l >= H->m) break;
        else if (r >= H->m) j = l;
        else {
            j = H->xs[l].count > H->xs[r].count ? l : r;
        }

        if (H->xs[i].count > H->xs[j].count) {
            hswap(H, i, j);
            i = j;
        }
        else break;
    }

    return true;
}


