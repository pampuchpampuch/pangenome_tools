import argparse
from pgtools import gff_parser, maf_parser, utils
import os
"""
summarises common contents of two MAFs. Uses bedtools to extract intersection.
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("maf1")
    parser.add_argument("maf2")
    parser.add_argument("out_dir")    
    args = parser.parse_args()
    # fasta = open(args.fasta,"w")
    # gff = open(args.gff, "r")
    # fasta = args.out_dir +"/" + args.gff.split("/")[-1][:-3]+"fa"
    # gff_res = open(os.path.join(args.out_dir, "all.gff"),"w")
    # genome_stats = {}
    mafs = [args.maf1, args.maf2]
    gff_paths = []
    for i in range(2):
        maf = maf_parser.parse_maf(mafs[i])
        maf.detect_soft_core()
        gff_path = os.path.join(args.out_dir, f"maf{i}.gff")
        maf.to_GFF(gff_path)
        gff_paths.append(gff_path)

    utils.bedtools_intersect_gff(gff_paths[0], gff_paths[1], out_dir=args.out_dir, out_gff="inter.gff")

    inter_gff  = os.path.join(args.out_dir, "inter.gff")
    gff_paths.append(inter_gff)
    row = mafs
    for gff in gff_paths:
        row.append(utils.get_gff_sum_aggregate_lens(gff))

    res_csv = open(os.path.join(args.out_dir, "summary.csv"))
    res_csv.write("maf1,maf2,maf1 pos n,maf2 pos n,common pos n\n")
    row_str = ",".join(row)
    res_csv.write(f"{row_str}\n")