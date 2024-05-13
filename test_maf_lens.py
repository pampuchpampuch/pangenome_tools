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
    
    args = parser.parse_args()

    maf = parse_maf(args.maf_in)

    simple_sum = 0
    for block in maf.synteny_blocks:
        for maf_seq in block.block_seqs:
            seq = "".join(maf_seq.seq.split("-"))
            len_seq = len(seq)
            simple_sum += len_seq
            len_maf = len(maf_seq)
            one_start, one_end = maf_seq.one_based_coords(maf_seq.start, maf_seq.end, maf_seq.chr_size, maf_seq.strand)
            len_one_based = one_end - one_start +1
            if not len_seq == len_maf == len_one_based:
                print(f"len_seq:{len_seq}")
                print(f"calculated",{len_maf})
                print(f"from one based: {len_one_based}")

        maf_seqs_gff_coords = {}
    
    maf_sum = 0
    maf_sequences = maf.get_sequences()
    for genome, seqs in maf_sequences.items():
        maf_seqs_gff_coords[genome] = sorted([seq.one_based_coords(seq.start, seq.end, seq.chr_size, seq.strand)
                                       for seq in seqs], key = lambda x: (x[0], x[1]))
    # maf.write_to_file(args.maf_out)
        
    assert simple_sum == maf_sum, f"{simple_sum}, {maf_sum}"

if __name__ == "__main__":
    main()