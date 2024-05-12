import argparse
import os
from pangenome_models.maf_parser import parse_maf

def fix_seq_len(maf_seq):
    seq = maf_seq.seq
    # seq_len = len("".join(seq.split("-")))
    seq = "".join(seq.split("-"))
    maf_seq.end = maf_seq.start + len(seq)
    # print("".join(seq.split("-")))
    # print(seq)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("maf_in",
                        help="path to first maf file")
    parser.add_argument("maf_out",
                        help="path to output maf file")
    
    args = parser.parse_args()

    maf = parse_maf(args.maf_in)

    for block in maf.synteny_blocks:
        for seq in block.block_seqs:
            print("stara",len(seq))
            fix_seq_len(seq)
            print("nowa",len(seq))

    maf.write_to_file(args.maf_out)

if __name__ == "__main__":
    main()