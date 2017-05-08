from __future__ import print_function
from __future__ import division

import pysam

import argparse

class Read(object):
    name = None
    max_quality = 0
    sequence = None

    def __init__(self, read, quality_threshold=40):
        self.name = read.query_name
        self.max_quality = read.mapping_quality
        if self.max_quality < quality_threshold:
            self.sequence = read.query_sequence

    def update(self, read, quality_threshold=40):
        assert self.name == read.query_name
        self.max_quality = max(self.max_quality, read.mapping_quality)
        if self.max_quality >= quality_threshold:
            self.sequence = None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_bam", help="Path to input bam")
    parser.add_argument("output_fastq", help="Path to output fastq")
    parser.add_argument("--quality_threshold", default=40, type=int, help="Read quality threshold (Default: %(default)s)")

    args = parser.parse_args()

    reads = pysam.AlignmentFile(args.input_bam, "rb")

    read_qual = {}

    for read in reads:
        if read.query_name in read_qual:
            read_qual[read.query_name].update(read)
        else:
            read_qual[read.query_name] = Read(read)
    
    with open(args.output_fastq, "w") as outfile:
        for name, read in read_qual.items():
            if read.max_quality < args.quality_threshold:
                print("@{}".format(read.name), file=outfile)
                print(read.sequence, file=outfile)
                print("+", file=outfile)
                print('"'*len(read.sequence), file=outfile)

