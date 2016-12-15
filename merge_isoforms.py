from __future__ import print_function
from __future__ import division

from Bio import SeqIO
import pandas as pd

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", help="Path to output table")
    parser.add_argument("input", nargs="+", help="List of input file names")

    args = parser.parse_args()

    dat_dict = {"group": [], "isoform": [], "length": []}
    for input in args.input:
        for record in SeqIO.parse(input, "fastq"):
            dat_dict["group"].append(record.id)
            dat_dict["isoform"].append(record.id)
            dat_dict["length"].append(len(record.seq))
  
    dat = pd.DataFrame.from_dict(dat_dict)
    dat.to_csv(args.output, sep="\t", index=False)
