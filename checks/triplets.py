from pgtools.utils import intersection_len, contains
from pgtools.gff_parser import parse_GFFs_dir
from pgtools.maf_parser import parse_maf, MAF
import copy 
from tqdm import tqdm
from pgtools import panaroo_parser

gff_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/gff"
maf_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/new_maf_no_refound.maf"
panaroo_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset/panaroo_out"
panaroo_obj = panaroo_parser.parse_panaroo_output(panaroo_dir, gff_dir, include_refound=True)
# simple_gffs = parse_GFFs_dir(gff_dir)
maf = parse_maf(maf_dir)


# #for maf from panaroo, precise coords should be found
# cds_strand = {}
# all_cds_coords = []
# for genome, gff in simple_gffs.simple_GFFs.items():
#     for scaff, cds in gff.scaffolds_coords.items():
#         for cds_ in cds:
#             all_cds_coords.append((cds_[0], cds_[1]))
#             cds_strand[(cds_[0], cds_[1])] = cds_[2]

# print(all_cds_coords)
# cds_overlaps = {}
# genome_overlaps = {}
# for genome, gff in simple_gffs.simple_GFFs.items():
#     gff_overlaps = gff.get_overlapping_annotations()
#     genome_overlaps[genome] = gff_overlaps
#     for scaff, overlaps_ in gff_overlaps.items():
#         cds_overlaps[f"{genome}.{scaff}"] = overlaps_

# print(cds_overlaps)

# for genome, overlaps_ in genome_overlaps.items():
#     print(genome)
#     print(overlaps_)
#     print("-"*40)        
# longest_overlaps = {}
# for genome, gff in simple_gffs.simple_GFFs.items():
#     longest_overlaps[genome] = gff.get_highest_overlap()

# print(longest_overlaps)


# for seq_coll in list(maf.seq_collections)[:3]:
#     print(seq_coll.id)
#     for seq in seq_coll.sequences:
#         print(seq.coord_system)
#     print("-"*30)
#     # break

# for seq_coll in list(panaroo_obj.seq_collections)[:3]:
#     print(seq_coll.id)
#     for seq in seq_coll.sequences:
#         print(seq.coord_system)
#     print("-"*30)
#     # break

panaroo_small = panaroo_parser.Panaroo(list(panaroo_obj.seq_collections)[:100], panaroo_obj.graph, panaroo_obj.soft_core_thresholds)

# print("panaroo")
# for seq_coll in panaroo_obj.seq_collections:
#     print(seq_coll.id)
#     break
#     # for seq in seq_coll.sequences:
#     #     # seq = seq
#     #     break

print("panaroo_small")
seq_panaroo = None
for seq_coll in panaroo_small.seq_collections:
    # print("panaroo_ids")
    # print(seq_coll.id)
    # print(seq_coll.cluster_name)
    for seq in seq_coll.sequences:
        seq_panaroo = seq
        break

genome = seq_panaroo.get_genome_name()
print("Genome", genome)
seq_name = seq_panaroo.seq_name

sequence_dict = panaroo_small.get_sequences_by_seq_name()
# print(sequence_dict)

panaroo_filtred = panaroo_small.filter_by_genome([genome])
print(panaroo_filtred.size())

clust_id_name_mapping = panaroo_filtred.get_clust_id_name_mapping()

# clust_adj_dict = panaroo_filtred.get_cluster_adjency_dict()
# print(clust_adj_dict)

# full_clust_adj_dict = panaroo_obj.get_cluster_adjency_dict()
# print("Full")
# len_counts = {}
# for k, v in full_clust_adj_dict.items():
#     # print(k, len(v))
#     # print(k,v)
#     if len(v) in len_counts:
#         len_counts[len(v)] += 1
#     else:
#         len_counts[len(v)] = 1

# print(len_counts)

# print("n_clusters", len(panaroo_obj.seq_collections))
        
"""
adj cluster dict for panaroo aln files does not exactly match gml from panaroo, but 
that is because Panaroo uses then multiple steps that modify the structure.
"""

panaroo_small.get_panaroo_triplets()
