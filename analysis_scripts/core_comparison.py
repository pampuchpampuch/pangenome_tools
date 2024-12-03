"""
Analysis consists of the following steps:
1. Save models after core detection as GFF files
2. Prepare aggreagted GFFs from panaroo input gffs
3. 

"""

import argparse
import os
import pandas as pd
from pgtools import maf_parser, gff_parser, utils

def clean_gff(gff, gff_res):
    """
    Writes gff without fasta part
    """
    # gff_res = open(os.path.join(out_dir, gff),"w")
    genome_name = gff.split("/")[-1][:-len(".gff")]
    with open(gff) as f:
        for line in f:
            if line.startswith("##FASTA"):
                # gff_res.close()
                return
            if line.startswith("#"):
                # gff_res.write(line)
                continue
            else:
                gff_res.write(f"{genome_name}.{line}")

def prepare_joined_gff(gff_dir, res_gff_dir):
    gff_res = open(res_gff_dir,"w")
    for f in os.listdir(gff_dir):
        f = os.path.join(gff_dir, f)
        if f.endswith(".gff"):
            clean_gff(f, gff_res)
    gff_res.close()

def get_only_core_seq(gff, res_gff):
    gff_res = open(res_gff,"w")
    with open(gff) as f:
        for line in f:
            if line.strip().endswith("True"):
                gff_res.write(line)
                # print(line)
    gff_res.close()

def get_cds_sum_lens(gff):
    cds_lens = {}
    gff = gff_parser.parse_joined_gff(gff)
    for cds in gff:
        cds_len = cds.end - cds.start + 1
        if cds.seq_name in cds_lens:
            if cds.annotation_id in cds_lens[cds.seq_name]:
                cds_lens[cds.seq_name][cds.annotation_id] += cds_len
            else:
                cds_lens[cds.seq_name][cds.annotation_id] = cds_len
        else:
            cds_lens[cds.seq_name] = {cds.annotation_id : cds_len}
    return cds_lens

def dict_to_ann_set(mapped_lens):
    res_cds = []
    for seq_name, mapped in mapped_lens.items():
        for ann_id, ann_len in mapped.items():
            res_cds.append(f"{seq_name}.{ann_id}")
    return set(res_cds)  
    
def filter_cds_lens(cds_lens, mapped_lens, thr=0.5):
    res_cds = []
    for seq_name in mapped_lens:
        full_lens = cds_lens[seq_name]
        mapped = mapped_lens[seq_name]
        for ann_id, ann_len in mapped.items():
            if ann_len >= thr * full_lens[ann_id]:
                res_cds.append(f"{seq_name}.{ann_id}")
    return set(res_cds)

def summarise_dataset(gff_all, gff_maf_prot, gff_core_prot, model_name): 

    gff_cds = gff_parser.parse_joined_gff(gff_all)
    all_cds_lens = {}
    for cds in gff_cds:
        if cds.seq_name in all_cds_lens:
            all_cds_lens[cds.seq_name][cds.annotation_id] = cds.end - cds.start + 1
        else:
            all_cds_lens[cds.seq_name] = {cds.annotation_id : cds.end - cds.start + 1}   
    core_lens = get_cds_sum_lens(gff_core_prot)
    prot_lens = get_cds_sum_lens(gff_maf_prot)
    
    signif_core = filter_cds_lens(all_cds_lens, core_lens)
    signif_prot = filter_cds_lens(all_cds_lens, prot_lens)
    all_core = dict_to_ann_set(core_lens)
    all_prot = dict_to_ann_set(prot_lens)
    source_all = dict_to_ann_set(all_cds_lens)

    return {model_name: {"cds_all": {"core": all_core, "all": all_prot}, "cds_significant": {"core": signif_core, "all": signif_prot}}}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("maf1")
    parser.add_argument("maf2")
    parser.add_argument("gffs")
    parser.add_argument("out_dir")
    parser.add_argument("--overlap_threshold", default=0.5, dest="overlap_threshold", type=float,
                        help="section of sequences that have to overlap to consider them significant")
    args = parser.parse_args()

    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # prepare gff with all proteins in the pangenome
    gff_joined_path = os.path.join(out_dir, "cds_all.gff")
    prepare_joined_gff(args.gffs, gff_joined_path)

    # save mafs as gff
    maf_paths = [args.maf1, args.maf2]
    for i in range(len(maf_paths)):
        maf_gff_path = os.path.join(out_dir, f"maf{i+1}.gff")
        maf_obj = maf_parser.parse_maf(maf_paths[i])
        maf_obj.detect_soft_core()
        maf_obj.to_GFF(maf_gff_path)

    # save maf gffs with only core sequences
    
    for i in range(len(maf_paths)):
        gff_path = os.path.join(out_dir, f"maf{i+1}.gff")
        core_gff_path = os.path.join(out_dir, f"maf{i+1}_core.gff")
        get_only_core_seq(gff_path, core_gff_path)

    # prepare intersecting gffs
    
    gff_prots = []
    gff_core = []
    for i in range(len(maf_paths)):
        core_inter_path = os.path.join(out_dir, f"maf{i+1}_core_prots.gff")
        core_gff_path = os.path.join(out_dir, f"maf{i+1}_core.gff")
        gff_core.append(core_inter_path)
        utils.bedtools_intersect_gff(gff_joined_path, core_gff_path, out_dir, f"maf{i+1}_core_prots.gff")

        inter_path = os.path.join(out_dir, f"maf{i+1}_prots.gff")
        gff_prots.append(inter_path)
        gff_path = os.path.join(out_dir, f"maf{i+1}.gff")
        utils.bedtools_intersect_gff(gff_joined_path, gff_path, out_dir, f"maf{i+1}_prots.gff")

    # now all needed gffs are ready. Below is statistics summary calculation
        
    gff_cds = gff_parser.parse_joined_gff(gff_joined_path)
    all_cds_lens = {}
    for cds in gff_cds:
        if cds.seq_name in all_cds_lens:
            all_cds_lens[cds.seq_name][cds.annotation_id] = cds.end - cds.start + 1
        else:
            all_cds_lens[cds.seq_name] = {cds.annotation_id : cds.end - cds.start + 1}   
    source_data_cds_n = dict_to_ann_set(all_cds_lens)
    summary_dict = {}
    for i in range(len(maf_paths)):
        summary_dict.update(summarise_dataset(gff_joined_path, gff_prots[i], gff_core[i], f"maf{i+1}"))

    cols_1 = ["maf 1 path", "maf 2 path", "maf 1 "]

    for model, summary in summary_dict.items():
        print(model)
        for cds_category, annots_n in summary.items():
            print(cds_category, len(annots_n))
            for annot_type, ann_ids in annots_n.items():
                print(annot_type, len(ann_ids))
    cols_1 = ["maf path", "model", "all annots n", "core annots n"]

    vals_1 = []
    for i in range(2):
        model = f"maf{i+1}"
        vals_1.append([maf_paths[i], model, len(summary_dict[model]["cds_significant"]["all"]), len(summary_dict[model]["cds_significant"]["core"])])

    df1 = pd.DataFrame(data=vals_1, columns=cols_1)
    df1.to_csv(os.path.join(out_dir, "models_protein_numbers.csv"), index=False)
    
    cols_2 = ["maf1", "maf2", "part of prots in maf1 that are in maf2", "part of prots in maf2 that are in maf1"]
    annots_1 = summary_dict["maf1"]["cds_significant"]["core"]
    annots_2 = summary_dict["maf2"]["cds_significant"]["core"]
    annots_common = annots_1.intersection(annots_2)
    vals_2 = [[maf_paths[0],
               maf_paths[1],
               len(annots_common) / len(annots_1),
               len(annots_common) / len(annots_2)]]
    
    df2 = pd.DataFrame(data=vals_2, columns = cols_2)
    df2.to_csv(os.path.join(out_dir, "core_proteins_stats.csv"), index=False)

if __name__ == "__main__":
    main()

