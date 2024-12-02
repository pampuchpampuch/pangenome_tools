import argparse
import os
from pgtools import gff_parser


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff_dir",
                        help="path to maf file")
    parser.add_argument("csv_out",
                        help="path to output gff file")
    
    args = parser.parse_args()
    # fasta = open(args.fasta,"w")
    # gff = open(args.gff, "r")
    # fasta = args.out_dir +"/" + args.gff.split("/")[-1][:-3]+"fa"
    gffs = gff_parser.parse_GFFs_dir(args.gff_dir, gff_simple=False)
    pg_gffs = gff_parser.Pangenome_Gffs_full(gffs)

    pg_gffs.scaff_lens_to_csv(args.csv_out)
if __name__ == "__main__":
    main()