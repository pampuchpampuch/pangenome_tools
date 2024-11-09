import argparse
import os

def clean_gff(gff, gff_res):
    """
    Writes gff without fasta part
    """
    # gff_res = open(os.path.join(out_dir, gff),"w")
    genome_name = gff.split("/")[-1][:-len(".gff")]
    with open(gff) as f:
        for line in f:
            if line.startswith("##FASTA"):
                # gff_res.close()
                return
            if line.startswith("#"):
                # gff_res.write(line)
                continue
            else:
                gff_res.write(f"{genome_name}.{line}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff",
                        help="path to dir with gff files with defined sequences with fasta headers at the end")
    parser.add_argument("out_dir",
                        help="path to new gff file")
    
    args = parser.parse_args()
    # fasta = open(args.fasta,"w")
    # gff = open(args.gff, "r")
    # fasta = args.out_dir +"/" + args.gff.split("/")[-1][:-3]+"fa"
    gff_res = open(os.path.join(args.out_dir, "all.gff"),"w")
    for f in os.listdir(args.gff):
        f = os.path.join(args.gff, f)
        if f.endswith(".gff"):
            clean_gff(f, gff_res)
    gff_res.close()
if __name__ == "__main__":
    main()