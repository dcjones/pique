
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include "version.h"

void print_help(FILE* fout)
{
    fprintf(fout,
"Usage: pique [option]... [file]...\n"
"Assemble short sequencing reads into contigs, take no prisoners.\n"
"Options:\n"
"  -h, --help           print this message\n"
"  -V, --version        display program version\n\n");

    /* TODO: An option to dump the de-bruijn graph rather that producing
     * contigs. */
}


int main(int argc, char* argv[])
{
    /* TODO:
     * 1. Read reads.
     * 2. Assemble.
     */


    return EXIT_SUCCESS;
}

