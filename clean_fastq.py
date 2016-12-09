from __future__ import print_function
from __future__ import division

import argparse

from Bio import SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Fastq file to remove duplicate records")
    parser.add_argument("output", help="Output fastq file")

    args = parser.parse_args()

    outfile = open(args.output, "w")

    ids = {}
    records = []
    for record in SeqIO.parse(args.input, "fastq"):
        if record.id not in ids:
            records.append(record)
        ids[record.id] = True
    SeqIO.write(records, args.output, "fastq")
        
