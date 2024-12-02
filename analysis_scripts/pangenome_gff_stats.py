import argparse
import os
from pgtools import gff_parser

def n50(cont_lens):
    n50 = None
    cont_lens.sort()
    # lens_rev = cont_lens[::-1]
    sum_lens = 0
    half_sum = sum(cont_lens) / 2
    i = 0
    while sum_lens < half_sum:
        sum_lens += cont_lens[i]
        n50 = cont_lens[i]
        i += 1
    return n50

def gff_stats(gff):
    """
    basic contig stats based on gff files from panaroo datasets
    """
    # gff_res = open(os.path.join(out_dir, gff),"w")\
    # genome_size = 0
    # contigs_n = 0
    # n50 = None
    genome_name = gff.split("/")[-1][:-len(".gff")]
    gff_info = gff_parser.parse_gff(gff)
    scaff_lens = [scaff.length for scaff in gff_info.scaffolds]

    return {genome_name: {"genome_size": sum(scaff_lens), "contigs_n": len(scaff_lens), "N50": n50(scaff_lens)}}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff",
                        help="path to dir with gff files with defined sequences with fasta headers at the end")
    parser.add_argument("csv_out",
                        help="csv summary")
    
    args = parser.parse_args()
    # fasta = open(args.fasta,"w")
    # gff = open(args.gff, "r")
    # fasta = args.out_dir +"/" + args.gff.split("/")[-1][:-3]+"fa"
    # gff_res = open(os.path.join(args.out_dir, "all.gff"),"w")
    genome_stats = {}

    for f in os.listdir(args.gff):
        f = os.path.join(args.gff, f)
        if f.endswith(".gff"):
            genome_stats.update(gff_stats(os.path.join(f)))
    
    cols = ["genome","genome_size","contigs_n","N50"]
    res_csv = open(args.csv_out, "w")
    res_csv.write(",".join(cols) + "\n")
    for k, v in genome_stats.items():
        stats = ",".join([str(v[i]) for i in cols[1:]])
        res_csv.write(f"{k},{stats}\n")
if __name__ == "__main__":
    main()