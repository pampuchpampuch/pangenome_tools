import argparse
import os
from pangenome_models import maf_parser
from pangenome_models import gff_parser

def parse_maf_comp(dir,file_name):
    """
    Parses mafComparator output
    """
    lines = open(file_name).readlines()
    csv_line =[dir]
    for i in range(len(lines)):
        if lines[i].strip().startswith("<homologyTests fileA"):
            csv_line.append(str(round(float(lines[i+2].split()[-1].split('"')[1]),4)))

    return csv_line

def maf_comp_summary(maf_comp_dir, out_csv):
    """
    Summarizes mafComparator outputs for files in subdirectories
    (dir name is saved in csv)
    
    Parameters:
        maf_comp_dir: str
            directory with directories containing mafComparator outputs
    """
    csv_file = open(out_csv,"w")
    csv_file.write("dataset,recall,precision\n")

    for dataset in os.listdir(maf_comp_dir):
        for f in os.listdir(os.path.join(maf_comp_dir, dataset)):
            #file names have to follow strict naming convention
            if not f.endswith("mafComp.xml"): continue
            f_path = os.path.join(maf_comp_dir, dataset, f)
            csv_line=",".join(parse_maf_comp(dataset, f_path))+"\n"
            csv_file.write(csv_line)

def parse_maf_stats(file_path, dataset):
    """
    Parses mafStats output
    """
    # file = os.path.join(dir,file_name)
    lines = open(file_path).readlines()
    maf_source = lines[0].split('.')[0].split('_')[0]
    maf_stats = ["char_n", "gap_n", "columns_n",
                 "blocks_n", "avg_area", "max_area",
                 "total_area", "avg_degree", "max_degree",
                 "avg_seq", "max_seq", "n_seq"]
    stats_dict = {"dataset": dataset, "maf_source": maf_source}

    stats = []
    for i in range(12,24):
        if i == 15: continue
        stats.append(lines[i].split(":")[1].split()[0])
    
    stats.append(lines[25].split()[0])
    for i in range(len(maf_stats)):
        stats_dict[maf_stats[i]] = float(stats[i])
    stats_dict["gap_pct"] = round(stats_dict["gap_n"] / (stats_dict["gap_n"] + stats_dict["char_n"]), 2)

    return stats_dict

def maf_stats_summary(maf_stats_dir, out_csv):
    """
    Summarizes mafStats ouput
    """
    csv_file = open(out_csv,"w")
    col_names = "dataset,maf_source,gap_pct,blocks_n,total_area,avg_degree,max_degree,avg_seq,max_seq,n_seq"
    csv_file.write(col_names + "\n")

    for dataset in os.listdir(maf_stats_dir):
        for f in os.listdir(os.path.join(maf_stats_dir, dataset)):
            #file names have to follow strict naming convention
            if not f.endswith("stats.txt"): continue
            f_path = os.path.join(maf_stats_dir, dataset, f)
            # dataset_path = os.path.join(f_path, dataset)
            # csv_line=",".join(parse_maf_stats(f_path, dataset))+"\n"
            maf_stats = parse_maf_stats(f_path, dataset)
            csv_line = ",".join([str(maf_stats[k]) for k in col_names.split(",")])
            csv_file.write(csv_line+"\n")
    
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
                    for coords in common_area_cds[1:-2]:
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

def maf_cds_summary(maf_dir, gffs_dir, out_csv):
    """
    Sumamrizes cds statictics for maf files is subdirectories.
    In csv file, both subdir and file will be present.
    """
    csv_file = open(out_csv, "w")
    csv_colnames = "dataset,maf_source,%_cds_scaffolds,%_cds_in_maf,%_maf_cds_in_cds\n"
    csv_file.write(csv_colnames)

    for dataset in os.listdir(maf_dir):
        for f in os.listdir(os.path.join(maf_dir, dataset)):
            if not f.endswith(".maf"): continue
            gff_dir = os.path.join(gffs_dir,dataset,"gff")
            maf_path = os.path.join(maf_dir,dataset,f)
            cds_stats = maf_cds_content(maf_path, gff_dir)
            maf_source = f.split("_")[0]
            stats = ",".join([str(stat) for stat in cds_stats.values()])
            csv_line = f"{dataset},{maf_source},{stats}\n"
            csv_file.write(csv_line)

