from pgtools.utils import intersection_len, contains
from pgtools.gff_parser import parse_GFFs_dir
from pgtools.maf_parser import parse_maf, MAF
import copy 
from tqdm import tqdm
from pgtools import panaroo_parser

gff_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/gff"
maf_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/new_maf_no_refound.maf"
panaroo_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset/panaroo_out"
# panaroo_obj = panaroo_parser.parse_panaroo_output(panaroo_dir, gff_dir, include_refound=False)
# simple_gffs = parse_GFFs_dir(gff_dir)
maf = parse_maf(maf_dir)
small_maf = MAF(list(maf.seq_collections)[:1])
# for block in small_maf.seq_collections:
#     print(block.to_MAF_block())

# small_maf = small_maf.convert_coords_system("gff")
small_maf.map_to_gff(gff_dir)
# print("="*100)
# for block in small_maf.seq_collections:
#     print(block.to_MAF_block())

for seq_coll in small_maf.seq_collections:
    for seq in seq_coll.sequences:
        for i in seq.mapped_annotations:
            print(i.annotation_id)