import argparse
import os
import maf_parser
import gff_parser
    
def gff_cds_pct(gff_path):
    """
    calculates %cds in scaffold in gff file

    Parameters:
    gff_path: str
        path to gff file

    Returns: float
        %of cds in scaffolds

    """
    scaffolds_len = 0
    cds_len = 0
    with open(gff_path) as f:
        next(f)
        genome_name = gff_path.split("/")[-1][:-4]
        for line in f:

            if line.startswith("##FASTA"):
                return round(cds_len/scaffolds_len,2), scaffolds_len
            
            if line.startswith("##sequence-region"):
                _, scaffold_name, _, scaff_len = line.split()
                scaffolds_len += int(scaff_len)

            else:
                scaffold_name, _, _, start, end, _, strand, _, gene_info, *_ =line.split()
                annotation_id = gene_info.split(";")[0][3:]

                cds_len += int(end) - int(start) + 1

    return round(cds_len/scaffolds_len,2), scaffolds_len

def maf_cds_content(maf, gff_dir):
    """
    Calculates % of cds area in maf((intersection of maf seqs and cds seqs) / sum of seqs in maf)

    Parameters:
        maf: path to maf file (seq naming convention: <genome_name>.<scaffold_name>)
        gff_dir: path to file with gff files corresponding to sequences in maf
    """
    maf = maf_parser.parse_maf(maf, store_seqs=False)
    gffs = [gff_parser.parse_gff(os.path.join(gff_dir,filename)) for filename in os.listdir(gff_dir)]
    
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
    return round(sum_maf_cds_lens/sum_maf_lens,4)