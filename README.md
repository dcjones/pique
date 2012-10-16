
# Pique

A stand-alone assembler based on the methods used for the
[Quip](https://github.com/dcjones/quip) compression algrothim.

Pique is a pretty standard de Bruijn graph assembler, but it uses probabalistic
data structures to achive extreme space efficiency, while remaining quite fast.

This is a research project, and at this point I do not recommend using it for
any real analysis, but if you're curious about how efficient assembly can be,
give it a try.

## Install

Download the latest tarball. Run `./configure && make install`. Or, clone the
git repository and run `autoreconf -i && ./configure && make install`.

## Usage

Simple running `pique reads.fastq > contigs.fasta` will build a de Bruijn graph
from the reads, then traverse it, outputing contigs in the fasta format to
`contigs.fasta`.

There are a number of options which you can read about with `pique --help`.

