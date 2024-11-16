import argparse
import os
from pgtools import gff_parser

def all_scaffs_to_gff(gffs, out_file):
    res_gff = open(out_file, "w")
    # gff_all = gff_parser.parse_GFFs_dir(gffs, gff_simple=False)
    for genome, gff in gffs.items():
        for scaff in gff.scaffolds:
            # cont_n += 1
            # cont_lens.append(scaff.length)
            seq_name = f"{scaff.genome}.{scaff.name}"
            start = 1
            end = scaff.length
            strand_sign = "+"
            res_gff.write(f"{seq_name}\tunidentified\tunidentified\t{start}\t{end}\t.\t{strand_sign}\t0\tINFO=whole_scaffold\n")

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
    
    # res_file = open(args.out_file, "w")
    # res_file.write(all_scaffs_to_gappy_block(args.gff_dir))
    # res_file.close()
    gffs = gff_parser.parse_GFFs_dir(args.gff_dir, gff_simple=False)
    all_scaffs_to_gff(gffs, args.out_file)

if __name__ == "__main__":
    main()