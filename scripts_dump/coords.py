import argparse
from pgtools.analysis import maf_cds_content
from pgtools import maf_parser

""" TO DO
1. Sprawdz czy długość są równe:
    - dł sekwencji w MSA z panaroo
    - dł wyliczana z gff
    - dł wyliczana z MAF
        - mafowe coords
        - coords po konwersji


"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("maf_in",
                        help="path to first maf file")
    parser.add_argument("gff_dir",
                        help="path to the gff file_dir")  
    
    args = parser.parse_args()

    maf = maf_parser.parse_maf(args.maf_in)

    # simple_sum = 0
    # for block in maf.synteny_blocks:
    #     for maf_seq in block.block_seqs:
    #         seq = "".join(maf_seq.seq.split("-"))
    #         len_seq = len(seq)
    #         simple_sum += len_seq
    #         len_maf = len(maf_seq)
    #         one_start, one_end = maf_seq.one_based_coords(maf_seq.start, maf_seq.end, maf_seq.chr_size, maf_seq.strand)
    #         len_one_based = one_end - one_start +1
    #         if not len_seq == len_maf == len_one_based:
    #             print(f"len_seq:{len_seq}")
    #             print(f"calculated",{len_maf})
    #             print(f"from one based: {len_one_based}")
    # print(f"len_seq:{simple_sum}")
    #     maf_seqs_gff_coords = {}
    
    # maf_sum = 0
    # maf_sequences = maf.get_sequences()
    # for genome, seqs in maf_sequences.items():
    #     maf_seqs_gff_coords[genome] = sorted([seq.one_based_coords(seq.start, seq.end, seq.chr_size, seq.strand)
    #                                    for seq in seqs], key = lambda x: (x[0], x[1]))
    # # maf.write_to_file(args.maf_out)
        
    # assert simple_sum == maf_sum, f"{simple_sum}, {maf_sum}"
    res = maf_cds_content(args.maf_in, args.gff_dir)
    print(res)
    
if __name__ == "__main__":
    main()