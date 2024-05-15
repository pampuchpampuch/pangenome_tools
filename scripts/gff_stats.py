import os
import argparse
from pgtools.analysis import gff_cds_pct
    
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("gff_dir",
                        help="path to dir with gff files used to generate panaroo output")
    parser.add_argument("stats_out",
                        help="path to dir to write stats to")
    args = parser.parse_args()

    gff_n = 0
    cds_pct = 0

    stats_file = open(args.stats_out,"w")
    stats_file.write("file,sum_scaff_sum_len,longest_contig,longest_contig_contr\n")

    scaff_dict = {}
    best_scaff_dict = {}

    for filename in os.listdir(args.gff_dir): 
        _cds_pct, scaff_len, best_scaff, max_scaff_pct = gff_cds_pct(os.path.join(args.gff_dir,filename))
        # gff_out = gff_cds_pct(os.path.join(args.gff_dir,filename))
        # print(gff_out)
        cds_pct += _cds_pct
        gff_n += 1
        # stats_file.write(f"{filename},{scaff_len}\n")
        scaff_dict[filename] = scaff_len
        best_scaff_dict[filename] = [best_scaff,max_scaff_pct]


    print(scaff_dict)
    
    gff_sorted = sorted(scaff_dict.items(), key=lambda kv: kv[1])

    for gff, scaff_len in gff_sorted:
        stats_file.write(f"{gff},{scaff_len},{best_scaff_dict[gff][0]},{best_scaff_dict[gff][1]}\n")

    print(cds_pct/gff_n)


if __name__ == "__main__":
    main()
            


