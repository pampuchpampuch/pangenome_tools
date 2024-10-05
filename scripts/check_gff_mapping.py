from pgtools.utils import intersection_len, contains
from pgtools.gff_parser import parse_GFFs_dir
from pgtools.maf_parser import parse_maf, MAF
import copy 
from tqdm import tqdm

gff_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/gff"
maf_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/new_maf_no_refound.maf"

simple_gffs = parse_GFFs_dir(gff_dir)
maf = parse_maf(maf_dir)


#for maf from panaroo, precise coords should be found
cds_strand = {}
all_cds_coords = []
for genome, gff in simple_gffs.simple_GFFs.items():
    for scaff, cds in gff.scaffolds_coords.items():
        for cds_ in cds:
            all_cds_coords.append((cds_[0], cds_[1]))
            cds_strand[(cds_[0], cds_[1])] = cds_[2]

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

maf_block = None
for seq_coll in maf.seq_collections:
    if 4 < len(seq_coll) < 6:
        maf_block = copy.deepcopy(seq_coll)
        break

print(len(maf_block))

small_maf = MAF([maf_block])
small_maf_maf = copy.deepcopy(small_maf)
small_maf.convert_coords_system("gff")
small_maf_gff = small_maf
# print(small_maf)
# why no sequences after conversion??

def gff_coords(start, end, chr_size, strand):
    """
    For now, converts to gff compatible coords
    """
    #TO DO: Converts zero based, half-open coords (used in maf) to one-based, fully closed
    # print(strand)
    if strand>0:
        one_start = start
        one_end = end - 1
    else:
        one_start = chr_size - end
        one_end = chr_size - start

    return(one_start, one_end)

print("MAF")
print(list(small_maf_maf.seq_collections)[0].to_MAF_block())
print("GFF")
print(list(small_maf_gff.seq_collections)[0].to_MAF_block())

block_maf_coords = list(small_maf_maf.seq_collections)[0]

ugly_coords = []
for seq in sorted(block_maf_coords.sequences, key = lambda x: x.seq_name):
    ugly_coords.append(gff_coords(seq.start, seq.end, seq.src_size, seq.strand))

# print("conversion ugly")
# print(ugly_coords)

block_gff_coords = list(small_maf_gff.seq_collections)[0]


# try and map each sequence with converted coords
aut_gff_coords = []
aut_gff_ver = []
strandness = []
aut_gff_ver_strand = []
for seq in sorted(block_gff_coords.sequences, key = lambda x: x.seq_name):
    seq_coords = (seq.start, seq.end)
    aut_gff_coords.append(seq_coords)
    print(seq_coords)
    print(seq_coords in all_cds_coords)
    aut_gff_ver.append(seq_coords in all_cds_coords)
    strandness.append(seq.strand)
    aut_gff_ver_strand.append(cds_strand.get(seq_coords, False))


ugly_gff_ver_strand = []
ugly_gff_ver = []
for seq_coords in ugly_coords:
    ugly_gff_ver.append(seq_coords in all_cds_coords)
    ugly_gff_ver_strand.append(cds_strand.get(seq_coords, False))

aut_maf_coords = []
for seq in sorted(block_maf_coords.sequences, key = lambda x: x.seq_name):
    seq_coords = (seq.start, seq.end)
    aut_maf_coords.append(seq_coords)

print(strandness)
print("gff ugly")
print(ugly_coords)
print(ugly_gff_ver)
print(ugly_gff_ver_strand)
print("gff aut")
print(aut_gff_coords)
print(aut_gff_ver)
print(aut_gff_ver_strand)
print("maf")
print(aut_maf_coords)

maf.convert_coords_system("gff")

# for seq_coll in maf.seq_collections:
    # print(seq_coll)

block = None

# for seq_coll in maf.seq_collections:
#     print(seq_coll)

maf_obj = parse_maf(maf_dir)
model = maf_obj
model.convert_coords_system("gff")
# old_model = model.convert_coords_system("gff")

err_seqs = []
print("Checking if all map")
for seq_coll in tqdm(model.seq_collections):
    for seq in seq_coll.sequences:
        # print("seq")
        # print(seq.coord_system)
        assert((seq.start, seq.end) in all_cds_coords)

print("All seqs")
# print(old_model.coord_system)


# for seq in block.sequences:
#     print(seq.coord_system)
