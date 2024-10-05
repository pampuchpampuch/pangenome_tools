from pgtools.utils import intersection_len, contains, strand_rep
from pgtools.gff_parser import parse_GFFs_dir
from pgtools.maf_parser import parse_maf, MAF
from pgtools.panaroo_parser import parse_panaroo_output
import copy 
from tqdm import tqdm

gff_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/gff"
maf_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/new_maf_no_refound.maf"
panaroo_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/panaroo_out"

# print(strand_rep(-1))
# print(strand_rep(1))

simple_gffs = parse_GFFs_dir(gff_dir)
maf = parse_maf(maf_dir)
maf.convert_coords_system("gff")

panaroo_obj = parse_panaroo_output(panaroo_dir, gff_dir)


#for maf from panaroo, precise coords should be found
# cds_strand = {}
all_cds_coords = {}
for genome, gff in simple_gffs.simple_GFFs.items():
    for scaff, cds in gff.scaffolds_coords.items():
        all_cds_coords[f"{genome}.{scaff}"] = set()
        for cds_ in cds:
            all_cds_coords[f"{genome}.{scaff}"].add(cds_)

# print(all_cds_coords)


def check_annots(model, all_cds_coords):
    for seq_coll in model.seq_collections:
        for seq in seq_coll.sequences:
            if all_cds_coords[seq.seq_name]:
                for annot in all_cds_coords[seq.seq_name]:
                    # print(annot)
                    if contains((seq.start, seq.end),(annot.start, annot.end), threshold=0.8):
                        # seq.annotation_ids.add(annot.annotation_id)
                        print(annot.annotation_id)

                print(seq.annotation_ids)

panaroo_old = copy.deepcopy(panaroo_obj)
# panaroo_obj.convert_coords_system("gff")
# check_annots(panaroo_obj, all_cds_coords)

maf = parse_maf(maf_dir)

panaroo_obj.map_to_gff(gff_dir, overlap_threshold=1)
# maf.map_to_gff(gff_dir, overlap_threshold=1)

"""
SOMETIMES MULTIPLE CDS MAPPED - TRY USING STRAND ALSO
"""