#!/usr/bin/env python

__version__ = '1.0'

from Bio import SeqIO
import os, sys

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser("Convert fastq to fasta")
    parser.add_argument("fasta_filename", help="input fastq (must end with .fastq or .fq)")
    args = parser.parse_args()

    input = args.fasta_filename

    fq2fa(input)

def fq2fa(input):
    try:
        assert input.lower().endswith('.fastq') or input.lower().endswith('.fq')
    except AssertionError:
        print >> sys.stderr, "Input {0} does not end with .fastq or .fq! Abort".format(input)
        sys.exit(-1)
    output = input.replace(".fastq", "").replace(".fq", "") + '.fasta'

    f = open(output, 'w')
    for r in SeqIO.parse(open(input),'fastq'):
        SeqIO.write(r, f, 'fasta')
    f.close()

    print >> sys.stderr, "Output written to", f.name

if __name__ == "__main__":
    main()
