from pgtools import maf_parser
# from pgtools.gff_parser import parse_GFFs_dir, Gff, parse_gff
from pgtools.utils import intersection_len
from pgtools.pangenome import Pangenome
from itertools import combinations
from copy import deepcopy
from pgtools.pangenome import Pangenome
import os
import argparse
from tqdm import tqdm
import numpy as np 

def trim_overlap(s1, s2):
    if intersection_len((s1.start, s1.end), (s2.start, s2.end)) == 0:
        # return s1, s2
        return
    else:
        inter_len = intersection_len((s1.start, s1.end), (s2.start, s2.end))
    # print(inter_len)

    shorter_seq, longer_seq = sorted([s1, s2], key = lambda x: len(x))
    seq_len = len(shorter_seq.seq.replace("-",""))

    old_short = deepcopy(shorter_seq)
    # there is overlap and we want to trim the shorter sequence: 
    if shorter_seq.start < longer_seq.start:
        # longer seq is behind shorter seq - we want to trimm the end of shorter seq
        shorter_seq.end -= (inter_len + 1)
        if shorter_seq.strand > 0:
            shorter_seq.seq = shorter_seq.seq.replace("-","")[:-inter_len - 1]
        else:
            shorter_seq.seq = shorter_seq.seq.replace("-","")[inter_len + 1 :]
    else:
        shorter_seq.start += (inter_len + 1)
        if shorter_seq.strand > 0:
            shorter_seq.seq = shorter_seq.seq.replace("-","")[inter_len + 1:]
        else:
            shorter_seq.seq = shorter_seq.seq.replace("-","")[: -inter_len - 1]

def clean_block(block):
    for s1, s2 in combinations(block.sequences,2):
        # print("i")
        if s1.seq_name != s2.seq_name:
            continue
        inter_len = intersection_len((s1.start, s1.end), (s2.start, s2.end))
        if inter_len > 0:
            # print(s1, s2)
            trim_overlap(s1,s2)

def trimm_overlaps_maf(pangenome_obj, return_trimmed_ids = False) -> Pangenome:
    old_coord_sys = pangenome_obj.coord_system
    pangenome_obj = pangenome_obj.convert_coords_system("gff")

    # all_sequences = []
    # for block in pangenome_obj.seq_collections:
    #     # print(block.id)
    #     for seq in block.sequences:
    #         all_sequences.append(seq)
    # overlapping_pairs = {}
    # for s1, s2 in combinations(all_sequences,2):
    #     if s1.seq_name != s2.seq_name or s1.strand != s2.strand:
    #         continue
    #     inter_len = intersection_len((s1.start, s1.end), (s2.start, s2.end))
    #     if inter_len > 0:
    #         # print(s1, s2)
    #         overlapping_pairs[(s1, s2)] = inter_len

    # print(len(overlapping_pairs))
    # clean within block overlaps
    all_sequences = []
    # blocks_subset = deepcopy(list(maf.seq_collections)[:500])
    for block in pangenome_obj.seq_collections:
        # print(block.id)
        for seq in block.sequences:
            all_sequences.append(seq)
            
    print("getting sequence pairs")
    all_pairs_idx = list(combinations(list(range(len(all_sequences))),2))

    # this implementation was really slow, only combinations within contigs needed, changing approach

    contig_sequences_dict = {}
    for seq_idx, seq in enumerate(all_sequences):
        if seq.seq_name in contig_sequences_dict:
            contig_sequences_dict[seq.seq_name].append(seq_idx)
        else:
            contig_sequences_dict[seq.seq_name] = [seq_idx]
    
    same_contig_pairs_idx = []
    for contig, seq_idxs in contig_sequences_dict.items():
        same_contig_pairs_idx += list(combinations(seq_idxs, 2))

    # #filtering, so only pairs from same contig are later checked
    # same_contig_pairs_idx = []
    # for s1_idx, s2_idx in all_pairs_idx:
    #     if all_sequences[s1_idx].seq_name == all_sequences[s2_idx].seq_name:
    #         same_contig_pairs_idx.append((s1_idx, s2_idx))

    print("Finding overlaps")
    overlapping_pairs = {}
    for s1_idx, s2_idx in tqdm(same_contig_pairs_idx):
        s1 = all_sequences[s1_idx]
        s2 = all_sequences[s2_idx]
        inter_len = 0
        if s1.seq_name != s2.seq_name:
            continue
        inter_len = intersection_len((s1.start, s1.end), (s2.start, s2.end))
        if inter_len > 0:
            # print(s1, s2)
            overlapping_pairs[(s1, s2)] = inter_len
    trimmed_blocks = []
    for seq_pair in overlapping_pairs:
        for s in seq_pair:
            trimmed_blocks.append(s.cluster_id) 

    # print("Checking and optionally trimming within block overlaps")
    # for block in tqdm(pangenome_obj.seq_collections):
    #     clean_block(block)
        
        
    # all_sequences = []
    # # blocks_subset = deepcopy(list(maf.seq_collections)[:500])
    # for block in pangenome_obj.seq_collections:
    #     # print(block.id)
    #     for seq in block.sequences:
    #         all_sequences.append(seq)
            
    # overlapping_pairs = {}
    # print("Finding overlaps after cleaning blocks")
    # for s1, s2 in tqdm(combinations(all_sequences,2)):
    #     inter_len = 0
    #     if s1.seq_name != s2.seq_name:
    #         continue
    #     inter_len = intersection_len((s1.start, s1.end), (s2.start, s2.end))
    #     if inter_len > 0:
    #         # print(s1, s2)
    #         overlapping_pairs[(s1, s2)] = inter_len
            
    sorted_pairs = sorted(list(overlapping_pairs), key=lambda x: overlapping_pairs[x], reverse=True)

    print("Trimming sequences overlaps")
    for s1, s2 in tqdm(sorted_pairs):
        trim_overlap(s1, s2)
        
    new_blocks_dict = {}
    for seq in all_sequences:
        if seq.cluster_id in new_blocks_dict:
            new_blocks_dict[seq.cluster_id].add(seq)
        else:
            new_blocks_dict[seq.cluster_id] = {seq}

    new_blocks = []
    for block in pangenome_obj.seq_collections:
        if block.id in new_blocks_dict:
            # print(block.id)
            block.sequences = new_blocks_dict[block.id]
        new_blocks.append(block)
    pangenome_obj_new = deepcopy(pangenome_obj)
    pangenome_obj_new.seq_collections = new_blocks
    pangenome_obj = pangenome_obj.convert_coords_system(old_coord_sys)
    pangenome_obj_new = pangenome_obj_new.convert_coords_system(old_coord_sys)

    if return_trimmed_ids:
        # trimmed_blocks = []
        # for seq_pair in overlapping_pairs:
        #     for s in seq_pair:
        #         trimmed_blocks.append(s.cluster_id)       
        return pangenome_obj_new, trimmed_blocks
    
    return pangenome_obj_new

def get_overlapping_pairs(maf):
    small_maf = maf.convert_coords_system("gff")
    all_sequences = []
    # blocks_subset = deepcopy(list(maf.seq_collections)[:500])
    for block in small_maf.seq_collections:
        # print(block.id)
        for seq in block.sequences:
            all_sequences.append(seq)
            
    overlapping_pairs = {}
    for s1, s2 in combinations(all_sequences,2):
        inter_len = 0
        if s1.seq_name != s2.seq_name:
            continue
        inter_len = intersection_len((s1.start, s1.end), (s2.start, s2.end))
        if inter_len > 0:
            # print(s1, s2)
            overlapping_pairs[(s1, s2)] = inter_len
    return overlapping_pairs

def main():
    parser = argparse.ArgumentParser(description="Calculates how frequently\
                                     gfa verticles are present in maf blocks (for two genomes)")
    parser.add_argument("maf",
                        help="path to the maf file")
    parser.add_argument("maf_out",
                        help="path to the output (trimmed) maf")
    args = parser.parse_args()
    maf = args.maf
    maf = maf_parser.parse_maf(maf, store_seqs=True)

    # maf = maf_parser.MAF(list(maf.seq_collections)[:500])
    maf = maf_parser.MAF(list(maf.seq_collections))


    # maf.to_MAF("first.maf")
    # small_maf = deepcopy(small_maf_0)
    maf_trimmed, trimmed_block_ids = trimm_overlaps_maf(maf, return_trimmed_ids=True)

    # overlap_1 = get_overlapping_pairs(small_maf_0)
    # print("over 0", len(overlap_1))

    # overlap_2 = get_overlapping_pairs(maf_trimmed)
    # print("over 1", len(overlap_2))

    # blocks_dict = maf_trimmed.get_seq_collections_dict()
    # print("n trimmed blocks", len(trimmed_block_ids))
    # trimmed_block = blocks_dict[trimmed_block_ids[0]]
    # print(trimmed_block.to_MAF_block(align = True))   
    # for id in trimmed_block_ids:
    #     print(blocks_dict[id].to_MAF_block(align = True))
    #     break

    # print(maf_trimmed.coord_system)
    blocks_to_align = []
    for block in tqdm(maf_trimmed.seq_collections):
        if block.id in trimmed_block_ids:
            blocks_to_align.append(block)
    print("Aligning blocks if needed")
    for block in tqdm(blocks_to_align):
        # if block.id in trimmed_block_ids:
        block.align_sequences()

    # for block in tqdm(maf_trimmed.seq_collections):
    #     if block.id in trimmed_block_ids:
    #         block.align_sequences()
            # seq_lens = np.array([len(seq.seq) for seq in block.sequences])
            # seqs_refound = [seq.refound for seq in self.sequences]
            # assert all(seq_lens == seq_lens[0]), f"{block.id},{seq_lens},"     

    maf_trimmed.to_MAF(args.maf_out)

if __name__ == "__main__":
    main()