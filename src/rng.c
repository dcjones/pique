
#include <stdlib.h>

#include "rng.h"
#include "misc.h"

/* This is an implementation of the complementary multiple with carry (CMWC)
 * pseudo random number generator. */

struct rng_t_
{
    uint32_t Q[4096];
    uint32_t c;
    uint32_t i;
};


rng_t* rng_alloc(uint32_t seed)
{
    const uint32_t phi = 0x9e3779b9;
    rng_t* rng = malloc_or_die(sizeof(rng_t));
    rng->c = 362436;
    rng->i = 4095;

    rng->Q[0] = seed;
    rng->Q[1] = seed + phi;
    rng->Q[2] = seed + 2 * phi;

    uint32_t i;
    for (i = 3; i < 4096; ++i) {
        rng->Q[i] = rng->Q[i - 3] ^ rng->Q[i - 2] ^ phi ^ i;
    }

    return rng;
}


void rng_free(rng_t* rng)
{
    free(rng);
}


uint32_t rng_get(rng_t* rng)
{
    uint64_t t, a = UINT64_C(18782);
    uint32_t x, r = 0xfffffffe;
    rng->i = (rng->i + 1) & 4095;
    t = a * rng->Q[rng->i] + rng->c;
    rng->c = t >> 32;
    x = t + rng->c;
    if (x < rng->c) {
        ++x;
        ++rng->c;
    }

    return (rng->Q[rng->i] = r - x);
}


double rng_get_double(rng_t* rng)
{
    return (double) rng_get(rng) / (double) UINT32_MAX;
}

