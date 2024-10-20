import argparse
from tqdm import tqdm
from itertools import combinations
from pgtools import gfa_parser
from pgtools import maf_parser
from pgtools.utils import intersection_len, contains

# def intersection_len(s1_coords, s2_coords):
#     """
#     seq2 is assumed to be from MAF and seq1 from gfa (checking how well
#     seq1 fits into seq2)
#     """
#     s1_start, s1_end = s1_coords
#     s2_start, s2_end = s2_coords
#     # if (seq1_coords[0] >= seq2_coords[0]) and (seq1_coords[1] <= seq2_coords[1])
#     if s1_end < s2_start or s1_start > s2_end:
#         return 0

#     if s1_start < s2_start:
#         inter_start = s2_start
#     else:
#         inter_start = s1_start    

#     if s1_end > s2_end:
#         inter_end = s2_end
#     else:
#         inter_end = s1_end
        
#     return inter_end - inter_start + 1

# def contains(s1_coords, s2_coords, threshold = 0.7):
#     """
#     Assumes s1 is from gfa and s2 from maf. The goal is to check if
#     gfa vertex is contained in maf block. Threshold is a fraction of 
#     gfa seq that is contained in maf block sequence
#     """
#     inter_len = intersection_len(s1_coords, s2_coords)
#     return inter_len >= threshold




def vertices_freq(maf_blocks, gfa_verticles, threshold = 0.7):
    """
    Calculates gfa vertices freq in maf blocks. Consideres only blocks and verticles
    that have both contigs in them. 
    """
    ### !!! works for only 2 contigs
    ### COULD BE IMPROVED BY SORTING COORDS
    ### ISSUE: in theory, one gfa V could map to multiple maf blocks?
    contained_lens = [[],[]]
    all_lens = []
    for strandness in maf_blocks.keys():
        for V in gfa_verticles[strandness]:
            in_block = False
            for block in maf_blocks[strandness]:
                gff_vert_len = V[0][1]-V[0][0]+1
                all_lens.append(gff_vert_len)
                inter_len_1 = intersection_len(V[0], block[0])
                inter_len_2 = intersection_len(V[1], block[1])
                contained_lens[0].append(inter_len_1)
                contained_lens[1].append(inter_len_2)
                """
                TU JEST ~BŁĄD   
                czasem gfa może być rozłożony na kilka bloków mafa - wtedy kilka pojedynczych fragmentów dzielimy
                kilka razy - efektywnie jeden gfa coveraga jest dzielony przez pełną długość kilka razy
                """
                # if (inter_len_1 >= threshold * gff_vert_len) and (inter_len_2 >= threshold * gff_vert_len):
                #     ### TODO how to average that
                #     # if our focus is particular block - averaging for two seqs is
                #     # a sensible approach. If our focus is parcitular contig - average for
                #     # contigs separately
                #     # contained_lens.append((inter_len_1 + inter_len_2) / 2)
                #     contained_lens[0].append(inter_len_1)
                #     contained_lens[1].append(inter_len_2)
                # print("v",V)
                # print("block", block)
            #     if contains(V[0], block[0]) and contains(V[1], block[1]):
            #         contained_lens.append(V[0][1]-V[0][0]+1)
            #         in_block = True
            # if not in_block:
            #     not_contained_lens.append(V[0][1]-V[0][0]+1)
    return contained_lens, all_lens
    # return contained_lens, not_contained_lens

def vertices_freq(maf_blocks, gfa_verticles, threshold = 0.7):
    """
    Calculates gfa vertices freq in maf blocks. Consideres only blocks and verticles
    that have both contigs in them. 
    """
    ### !!! works for only 2 contigs
    ### COULD BE IMPROVED BY SORTING COORDS
    ### ISSUE: in theory, one gfa V could map to multiple maf blocks?
    contained_lens = [[],[]]
    all_lens = []
    for strandness in maf_blocks.keys():
        for V in gfa_verticles[strandness]:
            # in_block = False
            for block in maf_blocks[strandness]:
                gff_vert_len = v
                """
                TU JEST ~BŁĄD   
                czasem gfa może być rozłożony na kilka bloków mafa - wtedy kilka pojedynczych fragmentów dzielimy
                kilka razy - efektywnie jeden gfa coveraga jest dzielony przez pełną długość kilka razy
                """
                # all_lens.append(gff_vert_len)
                inter_len_1 = intersection_len(V[0], block[0])
                inter_len_2 = intersection_len(V[1], block[1])
                contained_lens[0].append(inter_len_1)
                contained_lens[1].append(inter_len_2)
                # """
                # TU JEST ~BŁĄD   
                # czasem gfa może być rozłożony na kilka bloków mafa - wtedy kilka pojedynczych fragmentów dzielimy
                # kilka razy - efektywnie jeden gfa coveraga jest dzielony przez pełną długość kilka razy
                # """
            """
            rozwiązanie błędu
            jak już pozbiramy wszystkie możliwe overlapy - dodajemy do all lens
            MOŻE TU POWINNAŚ JUŻ UŚREDNIĆ???? i potem podzielić przez liczbę wierzchołków dla pary
            """
            all_lens.append(gff_vert_len)

                # if (inter_len_1 >= threshold * gff_vert_len) and (inter_len_2 >= threshold * gff_vert_len):
                #     ### TODO how to average that
                #     # if our focus is particular block - averaging for two seqs is
                #     # a sensible approach. If our focus is parcitular contig - average for
                #     # contigs separately
                #     # contained_lens.append((inter_len_1 + inter_len_2) / 2)
                #     contained_lens[0].append(inter_len_1)
                #     contained_lens[1].append(inter_len_2)
                # print("v",V)
                # print("block", block)
            #     if contains(V[0], block[0]) and contains(V[1], block[1]):
            #         contained_lens.append(V[0][1]-V[0][0]+1)
            #         in_block = True
            # if not in_block:
            #     not_contained_lens.append(V[0][1]-V[0][0]+1)
    return contained_lens, all_lens
    # return contained_lens, not_contained_lens

# def coocuring_contigs(gfa: gfa_parser.SimpleVertices, maf: maf_parser.MAF):
#     maf_contigs = set()
#     # for block in tqdm(list(maf.seq_collections)):
#     #     maf_contigs = maf_contigs.union(set(combinations(sorted(block.get_sequence_names()), 2)))
#     # print("maf_contigs", maf_contigs)
#     gfa_contigs = set()
#     # for V in tqdm(list(gfa.seq_collections)):
#     #    gfa_contigs = gfa_contigs.union(set((combinations(sorted(V.get_sequence_names()),2))))
#     print("dupa")
#     out = p_map(dupa, list(gfa.seq_collections))
#     print("dupa po")
#     print(out)
#     # inter_contigs = maf_contigs.intersection(gfa_contigs)
#     # returns coocuring contig pairs that exists in both files
#     return inter_contigs


# ======

def coocuring_contigs(gfa: gfa_parser.SimpleVertices, maf: maf_parser.MAF):
    maf_contigs = []
    for block in tqdm(list(maf.seq_collections)):
        maf_contigs += list(combinations(block.get_sequence_names(), 2))
    # print("maf_contigs", maf_contigs)
    maf_contigs = set([tuple(x) for x in map(sorted, maf_contigs)])
    gfa_contigs = []
    for V in tqdm(list(gfa.seq_collections)):
       gfa_contigs += list(combinations(V.get_sequence_names(),2))
    gfa_contigs = set([tuple(x) for x in map(sorted, gfa_contigs)])
    # out = p_map(list(gfa.seq_collections))
    # print(out)
    inter_contigs = maf_contigs.intersection(gfa_contigs)
    # returns coocuring contig pairs that exists in both files
    return inter_contigs

# ======

# def test_maf_filtering(maf):
#     maf = maf_parser.parse_maf(maf, store_seqs=False)
#     # print(maf.chr_names)
#     maf_contigs = set()
#     for block in tqdm(list(maf.synteny_blocks)):
#         maf_contigs = maf_contigs.union(set(combinations(sorted(block.get_contig_names()), 2)))
#     # maf_contigs = maf.chr_names
#     # print(maf_contigs)
#     contig_pair = list(maf_contigs)[3]
#     print(contig_pair)
#     maf_blocks = maf.get_filtered_blocks_by_strand(contig_pair)
#     print(maf_blocks)



def main():
    parser = argparse.ArgumentParser(description="Calculates how frequently\
                                     gfa verticles are present in maf blocks (for two genomes)")
    parser.add_argument("maf",
                        help="path to the maf file")
    parser.add_argument("gfa",
                        help="path to the gfa file")
    parser.add_argument("csv_out",
                        help="csv output file")
    args = parser.parse_args()
    gfa = args.gfa
    maf = args.maf
    gfa = gfa_parser.parse_gfa1(gfa)
    print("gfa parsed")
    maf = maf_parser.parse_maf(maf, store_seqs=False)
    common_contig_pairs = coocuring_contigs(gfa, maf)
    csv_out = open(args.csv_out, "w")
    csv_out.write("contig1,contig2\n")
    for c1,c2 in common_contig_pairs:
        csv_out.write(f"{c1},{c2}\n").flush()
    csv_out.close()
if __name__ == "__main__":
    main()