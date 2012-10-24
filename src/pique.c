
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include "dbg.h"
#include "fastq.h"
#include "misc.h"
#include "version.h"

void print_help(FILE* fout)
{
    fprintf(fout,
"Usage: pique [option]... [file]... > out.mm\n"
"Assemble short sequencing reads into contigs, take no prisoners.\n\n"
"By default, output is an adjacency matrix representation of the\n"
"De Bruijn graph in matrix market exchange format.\n\n"
"Options:\n"
"  -n                   maxmimum number of unique k-mers (larger numbers use\n"
"                       more memory but allow potentially more accurate assembly\n"
"                       (default: 100000000)\n"
"  -k                   k-mer size used by the de bruijn (default: 25)\n"
"  -t, --threads        number of threads to use (default: 1)\n"
"  -h, --help           print this message\n"
"  -V, --version        display program version\n\n");

    /* TODO: An option to dump the de-bruijn graph rather that producing
     * contigs. */
}


typedef struct pique_ctx_t_
{
    fastq_t* f;
    pthread_mutex_t* f_mutex;
    dbg_t* G;
} pique_ctx_t;


void* pique_thread(void* arg)
{
    pique_ctx_t* ctx = arg;
    seq_t* seq = seq_create();
    twobit_t* tb = twobit_alloc();
    rng_t* rng = rng_alloc(1234);
    bool r;

    while (true) {
        pthread_mutex_lock(ctx->f_mutex);
        r = fastq_read(ctx->f, seq);
        pthread_mutex_unlock(ctx->f_mutex);
        if (!r) break;

        /* TODO: remove sequences with Ns? */

        twobit_copy_str_n(tb, seq->seq.s, seq->seq.n);
        dbg_add_twobit_seq(ctx->G, rng, tb);
    }

    rng_free(rng);
    seq_free(seq);
    return NULL;
}


int main(int argc, char* argv[])
{
    int opt, opt_idx;

    /* Size of the graph structure. */
    size_t n = 100000000;

    /* K-mer size. */
    size_t k = 25;

    /* Number of threads. */
    size_t num_threads = 1;

    static struct option long_options[] =
    {
        {"threads", required_argument, NULL, 't'},
        {"verbose", no_argument,       NULL, 'v'},
        {"help",    no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };

    while (true) {
        opt = getopt_long(argc, argv, "n:k:t:vh", long_options, &opt_idx);
        if (opt == -1) break;

        switch (opt) {
            case 'n':
                n = strtoul(optarg, NULL, 10);
                break;

            case 'k':
                k = strtoul(optarg, NULL, 10);
                break;

            case 't':
                num_threads = strtoul(optarg, NULL, 10);
                break;

            case 'v':
                pique_verbose = true;
                break;

            case 'h':
                print_help(stdout);
                return EXIT_SUCCESS;

            case '?':
                return 1;

            default:
                abort();
        }
    }

    kmer_init();
    dbg_t* G = dbg_alloc(n, k);

    pthread_mutex_t f_mutex;
    pthread_mutex_init_or_die(&f_mutex, NULL);

    pthread_t* threads = malloc_or_die(num_threads * sizeof(pthread_t));

    pique_ctx_t ctx;
    ctx.G = G;
    ctx.f_mutex = &f_mutex;
    size_t i;

    if (optind >= argc) {
        ctx.f = fastq_create(stdin);
        for (i = 0; i < num_threads; ++i) {
            pthread_create(&threads[i], NULL, pique_thread, &ctx);
        }

        for (i = 0; i < num_threads; ++i) {
            pthread_join(threads[i], NULL);
        }

        fastq_free(ctx.f);
    } else {
        FILE* file;
        for (; optind < argc; ++optind) {
            file = fopen(argv[optind], "r");
            if (file == NULL) {
                fprintf(stderr, "Cannot open %s for reading.\n", argv[optind]);
                return EXIT_FAILURE;
            }
            ctx.f = fastq_create(file);

            for (i = 0; i < num_threads; ++i) {
                pthread_create(&threads[i], NULL, pique_thread, &ctx);
            }

            for (i = 0; i < num_threads; ++i) {
                pthread_join(threads[i], NULL);
            }

            fastq_free(ctx.f);
        }
    }

    dbg_dump(G, stdout, num_threads);

    pthread_mutex_destroy(&f_mutex);
    dbg_free(G);
    kmer_free();

    return EXIT_SUCCESS;
}

