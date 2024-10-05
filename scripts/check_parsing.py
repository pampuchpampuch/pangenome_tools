from pgtools import panaroo_parser
from pgtools.gff_parser import parse_gff
import os
from Bio.Seq import Seq

panaroo_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset/panaroo_out"
gffs_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset/gff"
panaroo_obj = panaroo_parser.parse_panaroo_output(panaroo_dir, gffs_dir, include_refound=True)

print(panaroo_obj.soft_core_thresholds)
for seq_clust in panaroo_obj.seq_collections:
    # if seq_clust.soft_core:
    print(seq_clust.cluster_name)
    print(seq_clust.soft_core)
    # print()
    # for seq in seq_clust.sequences:
    #     print(seq.in_soft_core)
    # print("-"*20)

