import argparse
import os

def clean_gff(gff, out_dir):
    """
    Writes gff without fasta part
    """
    gff_res = open(os.path.join(out_dir, gff),"w")
    genome_name = gff.split(".")[0]
    with open(gff) as f:
        for line in f:
            if line.startswith("##FASTA"):
                gff_res.close()
                return
            if line.startswith("#"):
                # gff_res.write(line)
                continue
            else:
                gff_res.write(f"{genome_name}.{line}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff",
                        help="path to gff file with defined sequences with fasta headers at the end")
    parser.add_argument("out_dir",
                        help="path to output fasta file")
    
    args = parser.parse_args()
    # fasta = open(args.fasta,"w")
    # gff = open(args.gff, "r")
    # fasta = args.out_dir +"/" + args.gff.split("/")[-1][:-3]+"fa"
    clean_gff(args.gff, args.out_dir)
if __name__ == "__main__":
    main()