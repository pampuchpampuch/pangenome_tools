from pgtools import panaroo_parser
from pgtools.gff_parser import parse_gff, parse_GFFs_dir, Pangenome_Gffs
from pgtools.gfa_parser import parse_gfa1
from pgtools.maf_parser import parse_maf
from pgtools.utils import intersection_len, contains
from pgtools.pangenome import Pangenome
import os
from Bio.Seq import Seq

panaroo_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset/panaroo_out"
gffs_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset/gff"
gfa_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/klebsiella.gfa"
maf_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/cactus.maf"
# panaroo_obj = panaroo_parser.parse_panaroo_output(panaroo_dir, gffs_dir, include_refound=True)
# gfa_obj = parse_gfa1(gfa_dir)
maf_obj = parse_maf(maf_dir)
model = maf_obj
model.detect_soft_core(threshold=0.9)

n_core = 0
n_all = 0
for seq_coll in model.seq_collections:
    n_all += 1
    if seq_coll.soft_core: n_core += 1
    # print(seq_coll.core)
    # print(len(seq_coll.sequences))

print("pct_core", n_core / n_all)
print(n_core, n_all)

"""
Mapping to gff:
change model coord sys to gff (then all coords are for + strand - no symmetrical invert necessery)
if seq coll is core:
    for each seq:
        try mapping to gff:
            - use threshold for containing
            - go through all cds in appropriate scaffold and choose one that fits best

"""

pangenome_GFFs = parse_GFFs_dir(gffs_dir, gff_simple=True)
print(pangenome_GFFs.simple_GFFs)
for genome, scaffolds_coords in pangenome_GFFs.simple_GFFs.items():
    print(genome)
    print(scaffolds_coords)
    print("="*40)

# convert coord system to gff
    
model.convert_coords_system("gff")

for seq_coll in model.seq_collections:
    for seq in seq_coll.sequences:
        print("seq")
        print(seq)

def map_to_gff(model: Pangenome, gffs: Pangenome_Gffs):
    for seq_coll in model.seq_collections:
        for seq in seq_coll.sequences:
            genome = seq.get_genome_name()
            scaff = seq.get_contig_name()
            seq_cds = pangenome_GFFs[genome][scaff]
            best_coverage = 0
            best_cds = None
            for cds in seq_cds:
                coverage = intersection_len((seq.start, seq.end), (cds[0], cds[1]))
