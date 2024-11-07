from pgtools import panaroo_parser
# from pgtools.gff_parser import parse_GFFs_dir, Gff, parse_gff
from pgtools.utils import intersection_len
from pgtools.pangenome import Pangenome
from itertools import combinations
from copy import deepcopy
from pgtools.pangenome import Pangenome
from pgtools.maf_parser import SyntenyBlock, MAFseq
import os
import argparse
from tqdm import tqdm
import numpy as np 

def maf_new_blocks(panaroo_obj, maf_out):
    edges = panaroo_obj.get_cluster_edges()
    clusters_dict = {clust.cluster_name : clust for clust in panaroo_obj.seq_collections}
    src_sizes = {}
    for seq_coll in panaroo_obj.seq_collections:
        for seq in seq_coll.sequences:
            if seq.seq_name not in src_sizes:
                src_sizes[seq.seq_name] = seq.src_size

    res_maf = open(maf_out, "w")
    new_blocks = {}

    clusts_dicts = {}
    print("summarising blocks")
    # edges = edges[:5]
    for c1, c2 in tqdm(edges):
        clusts = [clusters_dict[c1], clusters_dict[c2]]
        # clust_2 = clusters_dict[c2]
        common_seqs = clusts[0].get_sequence_names().intersection(clusts[1].get_sequence_names())
        # print(common_seqs)
        
        if not common_seqs:
            continue

        # for each sequence: [most expreme coords for one cluster, most extreme poses for the other cluster]
        clusts_seqs_dict = {s:[] for s in common_seqs}
        # print(clusts_seqs_dict)
        for i in range(len(clusts)):
        # c1_seqs = {s.seq_name:None for s in common_seqs}
        # c2_seqs = {s.seq_name:None for s in common_seqs}
            clust_seq = {}
            for s in clusts[i].sequences:
                if s.seq_name not in common_seqs:
                    continue
                if s.seq_name in clust_seq:
                    s_old = clust_seq[s.seq_name]
                    # if false is assign to s_old - there is more 
                    # than one seqs from one seq with different strandness
                    if not s_old:
                        continue
                    start,end,strand = s_old
                    if s.strand != strand:
                        clust_seq[s.seq_name] = False
                    if s_old[0] > s.start:
                        start = s.start
                    if s_old[1] < s.end:
                        end = s.end
                    clust_seq[s.seq_name] = (start, end, strand)
                else:
                    clust_seq[s.seq_name] = (s.start, s.end, s.strand)
            for s, coords in clust_seq.items():
                clusts_seqs_dict[s].append(coords)

        # print(clusts_seqs_dict)
        clusts_dicts[f"{c1}-{c2}"] = clusts_seqs_dict
        # PROBLEM - CO JEŚLI W JEDNYM KLASTRZE JEDNA SEKWENCJA MA MNIEJSZE POS A NA DRUGIM WIĘKSZE, A INNA SEKWENCJANA ODWÓR?? - RODZIAŁ NA KLASTRY WAŻNNY (NA KOŃCU MOŻNA PO PROSTU PODZIELIĆ SEKWENCJE NA DWIE GRUPY)
    print("clusts_dicts",len(clusts_dicts))     
                
        # for seq_name in common_seqs
            # co jeśli jest paralog - wtedy weź sekwencję bliżej tej drugiej
            # if one seq start is bigger than other end - this is the interval for complement
        # in new blocks - left for seqs with lower starts on the left, right for seqs with lower starts on the right

    # print(len(clusts_dicts))
    # print(clusts_dicts[0])

    new_blocks_all = []
    for ids, clust in tqdm(list(clusts_dicts.items())):
        # print(clust)
        new_blocks = [[],[]]
        for seq, coords in clust.items():
            c1,c2 = coords
            #strands have to be the same
            if c1[-1] != c2[-1]:
                continue
            # print(seq, c1, c2)

            coords.sort(key = lambda x: x[0])
            new_start = coords[0][1] + 1
            new_end = coords[1][0] - 1
            # if sequences from neighbouring blocks are one after another without a gap or are overlapping: skip
            if new_start + 1 >= new_end:
                continue
            # print(coords)
            # print(c1)
            new_seq = MAFseq(seq, new_start, new_end, c1[-1], src_sizes[seq], seq="-")
            if c1[0] > c2[0]:
                new_blocks[1].append(new_seq)
            else:
                new_blocks[0].append(new_seq)

        if len(new_blocks[0]) or len(new_blocks[1]):
            clust_ids = ids.split("-")
            for c_id in clust_ids:
                block = clusters_dict[c_id]
                res_maf.write(f"#ORIGINAL:{c_id}\n")
                res_maf.write(block.to_MAF_block())
                res_maf.flush()
        for B in new_blocks:
            # print(B)
            if len(B) > 1:
                syn_block = SyntenyBlock(B)
                new_blocks_all.append(syn_block)
                # print(SyntenyBlock(B).to_MAF_block())
                res_maf.write(f"#ids:{ids}\n")
                res_maf.write(syn_block.to_MAF_block())
                res_maf.flush()

        # CHECK IF FOR SURE THE NEW ONES DO NOT OVERLAP THE OLD ONES
        all_new_seqs = new_blocks[0] + new_blocks[1]
        prev_seqs = {}
        for seq in all_new_seqs:
            prev_seqs[seq.seq_name] = {"old":[], "new":[(seq.start, seq.end)]}
        for c_id in ids.split("-"):
            block = clusters_dict[c_id]
            for seq in block.sequences:
                if not seq.seq_name in prev_seqs:
                    continue
                prev_seqs[seq.seq_name]["old"].append((seq.start, seq.end))

        # print(prev_seqs)
        for seq_name, seqs in prev_seqs.items():
            for coords_1 in seqs["new"]:
                for coords_2 in seqs["old"]:
                    assert intersection_len(coords_1, coords_2) == 0, seq_name
    res_maf.close()
    return new_blocks_all
        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("panaroo_dir",
                        help="path to the panaroo output")
    parser.add_argument("gff_dir")
    parser.add_argument("maf_out")
    args = parser.parse_args()

    panaroo_obj = panaroo_parser.parse_panaroo_output(args.panaroo_dir, args.gff_dir)
    new_blocks = maf_new_blocks(panaroo_obj, args.maf_out)
    # edges = panaroo_obj.get_cluster_edges()
    # print(edges[:4])
    # clusters_dict = panaroo_obj.get_seq_collections_dict()
    

if __name__ == "__main__":
    main()