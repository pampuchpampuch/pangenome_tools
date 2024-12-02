import argparse
import os
from pgtools import maf_parser


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("maf",
                        help="path to maf file")
    parser.add_argument("out_gff",
                        help="path to output gff file")
    
    args = parser.parse_args()
    # fasta = open(args.fasta,"w")
    # gff = open(args.gff, "r")
    # fasta = args.out_dir +"/" + args.gff.split("/")[-1][:-3]+"fa"
    maf = maf_parser.parse_maf(args.maf)
    maf.detect_soft_core()
    maf.to_GFF(args.out_gff)
if __name__ == "__main__":
    main()