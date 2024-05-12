import argparse

def gff2fasta(gff, fasta):
    fasta = open(fasta,"w")
    fasta_part = False
    with open(gff) as f:
        for line in f:
            if fasta_part:
                fasta.write(line)
                continue
            if line.startswith(">"):
                fasta_part=True
                fasta.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff",
                        help="path to gff file with defined sequences with fasta headers at the end")
    parser.add_argument("out_dir",
                        help="path to output fasta file")
    
    args = parser.parse_args()
    # fasta = open(args.fasta,"w")
    # gff = open(args.gff, "r")
    fasta = args.out_dir +"/" + args.gff.split("/")[-1][:-3]+"fa"
    gff2fasta(args.gff, fasta)
if __name__ == "__main__":
    main()