from __future__ import print_function
from __future__ import division


import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")

    args = parser.parse_args()

    dat = pd.read_table(args.input)
    summary = pd.DataFrame(columns=["partition", "flnc_in", "transcript_out", "unincorporated_flnc"], index=dat.partition.unique())

    for part in dat.partition.unique():
        piece = dat.ix[(dat.partition == part) & (dat.is_fl == "Y")]
        flnc_in = piece.shape[0]
        transcript_out = len(piece.pbid.unique())
        unincorporated_flnc = piece["stat"].value_counts()["unmapped"]
        summary.iloc[part] = [part, flnc_in, transcript_out, unincorporated_flnc]

    summary.to_csv(args.output, sep="\t", index=False)
