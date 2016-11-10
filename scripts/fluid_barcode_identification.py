import argparse
from classify_barcode import Classifier, ChimeraDetectionOptions

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads_fn", default="ccs.fasta", help="Circular consensus reads fasta (Default: %(default)s)")
    parser.add_argument("--primer_fn_forward", default="custom_barcode_primers_forward.fa", help="Forward primer barcode fasta (Default: %(default)s)")
    parser.add_argument("--primer_fn_reverse", default="custom_barcode_primers_reverse.fa", help="Reverse primer barcode fasta (Default: %(default)s)")
    parser.add_argument("--ncpus", default = 12, type=int, help="Number of cpus (Default: %(default)s)")
    parser.add_argument("--nfl_fn_out", default="isoseq_nfl.fasta",  help="Non-full length read fasta output (Default: %(default)s)")
    parser.add_argument("--flnc_fn_out", default="isoseq_flnc.fasta",  help="Full-length non-chimaeric read fasta output (Default: %(default)s)")

    args = parser.parse_args()

    opt = ChimeraDetectionOptions(50, 30, 100, 50, 150, True)
    c = Classifier(opts=opt, reuse_dom=True, reads_fn=args.reads_fn, primer_fn_forward=args.primer_fn_forward, primer_fn_reverse=args.primer_fn_forward, cpus=args.ncpus, out_nfl_fn=args.nfl_fn_out, out_flnc_fn=args.flnc_fn_out)
    c.run()
