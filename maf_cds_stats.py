import argparse
import os
from pangenome_models import maf_parser
from pangenome_models import gff_parser

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff_dir",
                        help="path to the directory with gff files")
    parser.add_argument("maf",
                        help="path to the maf file")
    parser.add_argument("file_out",
                        help="path to the output file")
    
    args = parser.parse_args()

    maf = maf_parser.parse_maf(args.maf, store_seqs=False)

    gffs = [gff_parser.parse_gff(os.path.join(args.gff_dir,filename)) for filename in os.listdir(args.gff_dir)]
    
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

        try:
            cds_coords = cds_seqs_coords[contig]
        except:
            continue

        for maf in maf_coords:
            maf_start, maf_end = maf[0], maf[1]
            sum_maf_lens += maf_end - maf_start + 1
            first_cds_idx = None
            last_cds_idx = None

            for i, _cds_coords in enumerate(cds_coords):
                cds_start, cds_end = _cds_coords[0], _cds_coords[1]

                # maf start earlier
                # then cds has to start before maf ends
                # or
                # cds starts earlier
                # then maf has to start before cds ends          
                if ((cds_start >= maf_start and cds_start <= maf_end) or
                    (cds_start < maf_start and cds_end > maf_start)):
                    if not first_cds_idx:
                        first_cds_idx = last_cds_idx = i
                    else:
                        last_cds_idx = i

            if first_cds_idx:
                common_area_cds = cds_coords[first_cds_idx:last_cds_idx+1]
                first_cds_start = common_area_cds[0][0]
                last_cds_end = common_area_cds[-1][1]
                
                area_start = first_cds_start if first_cds_start < maf_start else maf_start
                area_end = last_cds_end if last_cds_end < maf_end else maf_end

                if len(common_area_cds) == 1:
                    sum_maf_cds_lens += area_end - area_start + 1
                
                else:
                    # add area from first cds:
                    sum_maf_cds_lens += common_area_cds[0][1] - area_start + 1
                    # add area from last cds:
                    sum_maf_cds_lens += area_end - common_area_cds[-1][0] + 1

                    for coords in common_area_cds[1:-2]:
                        sum_maf_cds_lens += coords[1] - coords[0]

            else:
                continue
    print(f"%CDS in maf {args.maf}: {round(sum_maf_cds_lens/sum_maf_lens,4)}")

if __name__ == "__main__":
    main()