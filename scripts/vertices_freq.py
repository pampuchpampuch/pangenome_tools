import argparse
from tqdm import tqdm
from itertools import combinations
from pgtools import gfa_parser
from pgtools import maf_parser

def intersection_len(s1_coords, s2_coords):
    """
    seq2 is assumed to be from MAF and seq1 from gfa (checking how well
    seq1 fits into seq2)
    """
    s1_start, s1_end = s1_coords
    s2_start, s2_end = s2_coords
    # if (seq1_coords[0] >= seq2_coords[0]) and (seq1_coords[1] <= seq2_coords[1])
    if s1_end < s2_start or s1_start > s2_end:
        return 0

    if s1_start < s2_start:
        inter_start = s2_start
    else:
        inter_start = s1_start    

    if s1_end > s2_end:
        inter_end = s2_end
    else:
        inter_end = s1_end
        
    return inter_end - inter_start + 1

def contains(s1_coords, s2_coords, threshold = 0.9):
    """
    Assumes s1 is from gfa and s2 from maf. The goal is to check if
    gfa vertex is contained in maf block. Threshold is a fraction of 
    gfa seq that is contained in maf block sequence
    """
    inter_len = intersection_len(s1_coords, s2_coords)
    return inter_len >= threshold

def vertices_freq(maf_blocks, gfa_verticles, threshold = 0.9):
    """
    Calculates gfa vertices freq in maf blocks
    """
    ### !!! works for only 2 contigs
    ### COULD BE IMPROVED BY SORTING COORDS
    contained_lens = []
    all_lens = []
    for strandness in maf_blocks.keys():
        for V in gfa_verticles[strandness]:
            in_block = False
            for block in maf_blocks[strandness]:
                all_lens.append(V[0][1]-V[0][0]+1)
                inter_len_1 = intersection_len(V[0], block[0])
                inter_len_2 = intersection_len(V[1], block[1])
                if (inter_len_1 >= threshold) and (inter_len_2 >= threshold):
                    ### TODO how to average that
                    # if our focus is particular block - averaging for two seqs is
                    # a sensible approach. If our focus is parcitular contig - average for
                    # contigs separately
                    contained_lens.append((inter_len_1 + inter_len_2) / 2)


                # print("v",V)
                # print("block", block)
            #     if contains(V[0], block[0]) and contains(V[1], block[1]):
            #         contained_lens.append(V[0][1]-V[0][0]+1)
            #         in_block = True
            # if not in_block:
            #     not_contained_lens.append(V[0][1]-V[0][0]+1)
    return contained_lens, all_lens
    # return contained_lens, not_contained_lens

def coocuring_contigs(gfa: gfa_parser.SimpleVertices, maf: maf_parser.MAF):
    maf_contigs = set()
    for block in tqdm(list(maf.seq_collections)[:100]):
        maf_contigs = maf_contigs.union(set(combinations(sorted(block.get_sequence_names()), 2)))
    # print("maf_contigs", maf_contigs)
    gfa_contigs = set()
    for V in tqdm(list(gfa.seq_collections)[:200]):
        gfa_contigs = gfa_contigs.union(set((combinations(sorted(V.get_sequence_names()),2))))
    inter_contigs = maf_contigs.intersection(gfa_contigs)
    return inter_contigs

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

def gfa_vs_maf(gfa, maf, csv_out):
    # print(csv_out)
    csv_file = open(csv_out, "w")
    csv_file.write("contig_1, contig_2, frequency\n")
    gfa = gfa_parser.parse_gfa1(gfa)
    maf = maf_parser.parse_maf(maf, store_seqs=False)
    common_contig_pairs = coocuring_contigs(gfa, maf)
    # maf_contigs = set(maf.chr_names)
    # gfa_contigs = set()
    # for V in gfa.vertices.values():
    #     for genome in V.seq_dict.keys():
    #         gfa_contigs.add(genome)
    # all_contig_names = maf_contigs.intersection(gfa_contigs)
    # print(all_contig_names)
    # contig_pairs = combinations(all_contig_names, 2)
    contig_pairs = common_contig_pairs

    for contig_pair in contig_pairs:
        print(contig_pair)
        gfa_verticles = gfa.get_filtered_vertices_by_strand(contig_pair, symmetrical_invert=True)
        maf_blocks = maf.get_filtered_vertices_by_strand(contig_pair, symmetrical_invert=True)
        contained_lens, all_lens = vertices_freq(maf_blocks, gfa_verticles)
        if not contained_lens:
            freq = 0
        else: 
            freq = sum(contained_lens)/sum(all_lens)
        print("==================")
        # print(not_contained_lens)
        # if contained_lens:
        csv_file.write(f"{contig_pair[0]}, {contig_pair[1]}, {freq}\n")

    # # for contig_pair in contig_pairs:
    # gfa_verticles = gfa.get_filtered_vertices_by_strand(contig_pair)
    # maf_blocks = maf.get_filtered_blocks_by_strand(contig_pair)

    # contained_lens, not_contained_lens = vertices_freq(maf_blocks, gfa_verticles)
    # if not contained_lens:
    #     freq = 0
    # else: 
    #     freq = sum(contained_lens)/sum(contained_lens.extend(not_contained_lens))
    # # print(not_contained_lens)
    # if contained_lens:
    #     print(f"{contig_pair[0]}, {contig_pair[1]}, {freq}\n")

def test():
    maf_mock = [((100,500),(110,550)), ((1000,2000),(1010,2090))]

    maf_mock_coords =  {"++": maf_mock,
                        "+-": [],
                        "--": [],
                        "-+": []}

    gfa_mock_coords =  {"++": [((34,67),(145,490)), ((110,300),(120,400))],
                        "+-": [],
                        "--": [],
                        "-+": []}
    # for contig_pair in contig_pairs:
    gfa_verticles = gfa_mock_coords
    #gfa.get_filtered_vertices_by_strand(contig_pair)
    maf_blocks = maf_mock_coords
    #maf.get_filtered_blocks_by_strand(contig_pair)

    contained_lens, not_contained_lens = vertices_freq(maf_blocks, gfa_verticles)
    print(contained_lens)
    print(not_contained_lens)
    if not contained_lens:
        freq = 0
    else: 
        freq = sum(contained_lens)/sum(contained_lens + not_contained_lens)
    # print(not_contained_lens)
    print(freq)

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
    # test()
    # test_maf_filtering(args.maf)
    gfa_vs_maf(args.gfa, args.maf, args.csv_out)

if __name__ == "__main__":
    main()