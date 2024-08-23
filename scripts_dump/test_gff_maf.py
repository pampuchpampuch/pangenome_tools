import argparse
import os
from pgtools import maf_parser
from pgtools import gff_parser
from pgtools import analysis

def maf_cds_content(maf, gff_dir):
    """
    Calculates % of cds area in maf((intersection of maf seqs and cds seqs) / sum of seqs in maf)

    Parameters:
        maf: str
            path to maf file (seq naming convention: <genome_name>.<scaffold_name>)
        gff_dir: str
            path to file with gff files corresponding to sequences in maf
    
    Returns: Dict[str, int]
            dictionary with basic length info for maf and cds sequences
    """
    maf = maf_parser.parse_maf(maf, store_seqs=False)
    gffs = [gff_parser.parse_gff(os.path.join(gff_dir,filename)) for filename in os.listdir(gff_dir)]

    sum_scaff_lens = 0
    sum_cds_lens = 0

    for gff in gffs:
        sum_scaff_lens += sum([scaff.length for scaff in gff.scaffolds])
        sum_cds_lens += sum([len(cds) for cds in gff.cds])

    maf_sequences = maf.get_sequences()

    cds_sequences = {}

    for gff in gffs:
        cds_sequences.update(gff.get_sequences())
    
    cds_seqs_coords = {}
    for genome, seqs in cds_sequences.items():
        cds_seqs_coords[genome] = sorted([(seq.start, seq.end) for seq in seqs], key = lambda x: (x[0], x[1]))

    maf_seqs_gff_coords = {}
    for genome, seqs in maf_sequences.items():
        maf_seqs_gff_coords[genome] = sorted([seq.one_based_coords(seq.start, seq.end, seq.chr_size, seq.strand)
                                       for seq in seqs], key = lambda x: (x[0], x[1]))

    # find, how long are the fragments in maf sequences that are in the cds fragments

    sum_maf_lens = 0
    sum_maf_cds_lens = 0

    for contig in maf_seqs_gff_coords:
        maf_coords = maf_seqs_gff_coords[contig]
        for maf in maf_coords:
            maf_start, maf_end = maf[0], maf[1]
            sum_maf_lens += maf_end - maf_start + 1

        try:
            cds_coords = cds_seqs_coords[contig]
        except:
            continue

        for maf in maf_coords:
            maf_start, maf_end = maf[0], maf[1]
            # print("="*30)
            # print(f"maf:{maf_start}-{maf_end}")

            first_cds_idx = None
            last_cds_idx = None

            for i, _cds_coords in enumerate(cds_coords):
                cds_start, cds_end = _cds_coords[0], _cds_coords[1]
                # print(f"cds:{cds_start}-{cds_end}")


                # maf start earlier
                # then cds has to start before maf ends
                # or
                # cds starts earlier
                # then maf has to start before cds ends          
                if ((cds_start >= maf_start and cds_start < maf_end) or
                    (cds_start < maf_start and cds_end > maf_start)):
                    if not first_cds_idx:
                        first_cds_idx = i
                        last_cds_idx = i
                    else:
                        last_cds_idx = i

            if first_cds_idx != None:
                common_area_cds = cds_coords[first_cds_idx:last_cds_idx+1]
                # print(common_area_cds)
                first_cds_start = common_area_cds[0][0]
                last_cds_end = common_area_cds[-1][1]
                area_start = first_cds_start if first_cds_start > maf_start else maf_start
                area_end = last_cds_end if last_cds_end < maf_end else maf_end
                # print(area_start,area_end)
                if len(common_area_cds) == 1:
                    # print(sum_maf_cds_lens)
                    sum_maf_cds_lens += area_end - area_start + 1
                else:
                    # add area from first cds
                    sum_maf_cds_lens += common_area_cds[0][1] - area_start + 1
                    # add area from last cds
                    sum_maf_cds_lens += area_end - common_area_cds[-1][0] + 1
                    # add area for middle cds
                    for coords in common_area_cds[1:-1]:
                        sum_maf_cds_lens += coords[1] - coords[0] + 1
            else:
                continue

    # print("maf_lens",sum_maf_lens)
    # print("cds_lens",sum_maf_cds_lens)
    results = {"%_cds_scaffolds": sum_cds_lens/sum_scaff_lens,
              "%_cds_in_maf": sum_maf_cds_lens/sum_maf_lens,
              "%_maf_cds_in_cds": sum_maf_cds_lens/sum_cds_lens,
              }
    results = {k:round(v,4) for k,v in results.items()}

    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("maf_in",
                        help="path to first maf file")
    parser.add_argument("gff_dir",
                        help="path to the gff file_dir")  
    
    args = parser.parse_args()

    maf = maf_parser.parse_maf(args.maf_in)

    simple_sum = 0
    for block in maf.synteny_blocks:
        for maf_seq in block.block_seqs:
            seq = "".join(maf_seq.seq.split("-"))
            len_seq = len(seq)
            simple_sum += len_seq
            len_maf = len(maf_seq)
            one_start, one_end = maf_seq.one_based_coords(maf_seq.start, maf_seq.end, maf_seq.chr_size, maf_seq.strand)
            len_one_based = one_end - one_start +1
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
    print(f"len_seq:{simple_sum}")
    print(res)
    
if __name__ == "__main__":
    main()