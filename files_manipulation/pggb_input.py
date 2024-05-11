#dir with fasta files to pggb input - one fasta with all sequences with headers prepended with "[filename]."

import argparse
import os
from Bio import SeqIO

def concatenate_fastas(fasta_dir, out_dir, out_filename):

    if not os.path.exists(out_dir): os.makedirs(out_dir)

    out_file_path = os.path.join(out_dir, out_filename)
    out_file = open(out_file_path,"w")

    for filename in os.listdir(fasta_dir):
        if filename.endswith(".fa"):
            in_fasta = os.path.join(fasta_dir,filename)
            for record in SeqIO.parse(in_fasta, "fasta"):
                record.id  = filename[:-len(".fa")] + "." + record.id
                record.description = ""
                out_file.write(record.format("fasta"))

    out_file.close()    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_dir",
                        help="path to dir with fasta files")
    parser.add_argument("outdir",
                        help="output dir")
    parser.add_argument("out_file",
                        help="fasta to write fastas content to")
    
    args = parser.parse_args()
    
    concatenate_fastas(args.fasta_dir, args.outdir, args.out_file)

if __name__ == "__main__":
    main()