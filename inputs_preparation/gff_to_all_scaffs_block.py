import argparse
import os
from pgtools import gff_parser

def all_scaffs_to_gappy_block(gff_dir):
    pg_gffs = gff_parser.parse_GFFs_dir(gff_dir, gff_simple=False)
    scaffs = {}
    for genome, gff in pg_gffs.items():
        for scaff in gff.scaffolds:
            scaffs[scaff.genome + "." + scaff.name] = {"len": scaff.length, "first_pos": scaff.seq[0]}
    # return scaffs
    lines = "a\n"
    pre_gaps = 0
    post_gaps = len(scaffs) - 1
    for scaff_name, params in scaffs.items():
        line_seq = "-"*pre_gaps + params["first_pos"] + "-"*post_gaps
        pre_gaps += 1
        post_gaps -= 1
        strand_sign = "+"
        scaff_len = params["len"]
        lines += f"s\t{scaff_name}\t{0}\t{1}\t{strand_sign}\t{scaff_len}\t{line_seq}\n"
    return lines

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff_dir",
                        help="path to dir with gff files with defined sequences with fasta headers at the end")
    parser.add_argument("out_file",
                        help="path to out file")
    
    args = parser.parse_args()
    # fasta = open(args.fasta,"w")
    # gff = open(args.gff, "r")
    # fasta = args.out_dir +"/" + args.gff.split("/")[-1][:-3]+"fa"
    
    res_file = open(args.out_file, "w")
    res_file.write(all_scaffs_to_gappy_block(args.gff_dir))
    res_file.close()
if __name__ == "__main__":
    main()