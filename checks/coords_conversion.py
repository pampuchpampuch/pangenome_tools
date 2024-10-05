from pgtools.utils import intersection_len, contains
from pgtools.gff_parser import parse_GFFs_dir
from pgtools.maf_parser import parse_maf, MAF
import copy 
from tqdm import tqdm
from pgtools import panaroo_parser

gff_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/gff"
maf_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/new_maf_no_refound.maf"
panaroo_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset/panaroo_out"
panaroo_obj = panaroo_parser.parse_panaroo_output(panaroo_dir, gff_dir, include_refound=False)
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


for seq_coll in list(maf.seq_collections)[:3]:
    print(seq_coll.id)
    for seq in seq_coll.sequences:
        print(seq.coord_system)
    print("-"*30)
    # break

for seq_coll in list(panaroo_obj.seq_collections)[:3]:
    print(seq_coll.id)
    for seq in seq_coll.sequences:
        print(seq.coord_system)
    print("-"*30)
    # break

maf_gff = maf.convert_coords_system("gff")
maf_block = None
for seq_coll in list(maf_gff.seq_collections)[:3]:
    print(seq_coll.id)
    # if 4 < len(seq_coll) < 6:
    #     maf_block = copy.deepcopy(seq_coll)
    #     break
    for seq in seq_coll.sequences:
        # print(seq)
        print(seq.coord_system)
    print("-"*30)
    # break

maf_block = None
for seq_coll in maf.seq_collections:
    if 4 < len(seq_coll) < 6:
        maf_block = copy.deepcopy(seq_coll)
        break

# maf_block = None
# for seq_coll in maf.seq_collections:
#     # if 4 < len(seq_coll) < 6:
#     #     maf_block = copy.deepcopy(seq_coll)
#     #     break
#     for seq in seq_coll.sequences:
#         print(seq.coord_system)
#     break

# print(len(maf_block))

small_maf = MAF([maf_block])
small_maf_gff = small_maf.convert_coords_system("gff")

# small_maf_maf = copy.deepcopy(small_maf)
# small_maf_gff = small_maf
# # print(small_maf)
# # why no sequences after conversion??



print("MAF")
print(list(small_maf.seq_collections)[0].to_MAF_block())
print("GFF")
print(list(small_maf_gff.seq_collections)[0].to_MAF_block())



# print("conversion ugly")
# print(ugly_coords)

# block_gff_coords = list(small_maf_gff.seq_collections)[0]


# # try and map each sequence with converted coords


# maf.convert_coords_system("gff")

# # for seq_coll in maf.seq_collections:
#     # print(seq_coll)

# block = None

# # for seq_coll in maf.seq_collections:
# #     print(seq_coll)

# maf_obj = parse_maf(maf_dir)
# model = maf_obj
# model.convert_coords_system("gff")
# # old_model = model.convert_coords_system("gff")



# for seq in block.sequences:
#     print(seq.coord_system)
