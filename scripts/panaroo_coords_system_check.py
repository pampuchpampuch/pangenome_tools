from pgtools import panaroo_parser
from pgtools.gff_parser import parse_gff
import os
from Bio.Seq import Seq

panaroo_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset/panaroo_out"
gffs_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset/gff"
panaroo_obj = panaroo_parser.parse_panaroo_output(panaroo_dir, gffs_dir, include_refound=True)

# for clust in panaroo_obj.seq_collections:
#     print(clust.cluster_name)

# print(panaroo_obj.coord_system)
# print(panaroo_obj.graph)
# coll = panaroo_obj.seq_collections.pop()
# for seq in coll.sequences:
    # print(seq.start)
test_seqs = {"generic": {"+": None, "-": None}, "refound": {"+": None, "-": None}}

for seq_col in panaroo_obj.seq_collections:
    for gene in seq_col.sequences:
        if gene.refound:
            if gene.strand > 0:
                test_seqs["refound"]["+"] = gene
            else:
                test_seqs["refound"]["-"] = gene
        else:
            if gene.strand > 0:
                test_seqs["generic"]["+"] = gene
            else:
                test_seqs["generic"]["-"] = gene

print(test_seqs)            


gene_data = panaroo_parser.parse_gene_data(os.path.join(panaroo_dir, "gene_data.csv"))
csv_dict = {}
for csv_gene in gene_data.values():
    csv_dict[csv_gene.annotation_id] = csv_gene

for gene_type in test_seqs:
    for strand in ["+", "-"]:
        # print(test_seqs[gene_type][strand].seq_name)
        gene = test_seqs[gene_type][strand]
        gene_gff_file = os.path.join(gffs_dir, gene.get_genome_name()+".gff")
        gff_obj = parse_gff(gene_gff_file, store_sequences=True)
        # compare: cds seq, seq from scaffold splicing, seq from alignment
        # cds_sequences = gff_obj.get_sequences()[gene.seq_name]
        # cds_seq = gff_obj.get_cds_sequence(gene.annotation_id)
        scaff = gff_obj.get_scaffolds_dict()[gene.seq_name]
        csv_gene_seq = csv_dict.get(gene.annotation_id, "").seq
        csv_gene_coords = csv_dict.get(gene.annotation_id, "").refound_coords
        gene.seq = gene.seq.replace("-", "")
        # print(gene.refound)
        # print("scaff_len", len(scaff.seq))
        # print("csv coords", csv_gene_coords)
        # # print("scaff", scaff.seq[:30])
        # print(gene.annotation_id, gene.start, gene.end, gene.strand)
        # print("gff")
        # print("panaroo_aln", gene.seq)
        # print("csv", csv_gene_seq)
        # scaff_gene_seq = scaff.seq[gene.start:gene.end]
        scaff_gene_seq = Seq(scaff.seq)
        if gene.strand < 0:
            scaff_gene_seq = scaff_gene_seq.reverse_complement()
        scaff_gene_seq = scaff_gene_seq[gene.start:gene.end]
        # print("scaff_spliced", scaff_gene_seq)
        print(gene.seq.lower() == csv_gene_seq.lower() == scaff_gene_seq.lower())
        if not (gene.seq.lower() == csv_gene_seq.lower() == scaff_gene_seq.lower()):
            # incorrect_n += 1
            cds = gff_obj.get_cds_dict().get(gene.annotation_id, "NOT IN GFF")
            print(gene.annotation_id)
            print("refound", gene.refound)
            if not gene.refound:
                print("cds_coords", cds.start, cds.end)
                print("cds strand", cds.strand)
            print("gene_coords", gene.start, gene.end)
            print("strand", gene.strand)
            print("panaroo_aln\n", gene.seq.lower())
            print("csv\n", csv_gene_seq.lower())
            print("scaff\n", scaff_gene_seq.lower())
            print("="*30)
        # print("="*30)

all_refound = []
for seq_col in panaroo_obj.seq_collections:
    for gene in seq_col.sequences:
        if gene.refound:
            all_refound.append(gene)

for gene in all_refound:
    gene_gff_file = os.path.join(gffs_dir, gene.get_genome_name()+".gff")
    gff_obj = parse_gff(gene_gff_file, store_sequences=True)
    # compare: cds seq, seq from scaffold splicing, seq from alignment
    # cds_sequences = gff_obj.get_sequences()[gene.seq_name]
    # cds_seq = gff_obj.get_cds_sequence(gene.annotation_id)
    scaff = gff_obj.get_scaffolds_dict()[gene.seq_name]
    csv_gene_seq = csv_dict.get(gene.annotation_id, "").seq
    csv_gene_coords = csv_dict.get(gene.annotation_id, "").refound_coords
    gene.seq = gene.seq.replace("-", "")
    scaff_gene_seq = Seq(scaff.seq)
    if gene.strand < 0:
        scaff_gene_seq = scaff_gene_seq.reverse_complement()
    scaff_gene_seq = scaff_gene_seq[gene.start:gene.end]
    if not (gene.seq.lower() == csv_gene_seq.lower() == scaff_gene_seq.lower()):
        print(gene.refound)
        print("scaff_len", len(scaff.seq))
        print("csv coords", csv_gene_coords)
        # print("scaff", scaff.seq[:30])
        print(gene.annotation_id, gene.start, gene.end, gene.strand)
        print("panaroo_aln", gene.seq)
        print("csv", csv_gene_seq)
        # scaff_gene_seq = scaff.seq[gene.start:gene.end]
        print("scaff_spliced", scaff_gene_seq)
        print(gene.seq.lower() == csv_gene_seq.lower() == scaff_gene_seq.lower())
        print(len(gene.seq.lower()) == len(csv_gene_seq.lower()) == len(scaff_gene_seq.lower()))
        print("="*30)  

# gff_example = "5151_2_6.gff"

# example_gff = parse_gff(os.path.join(gffs_dir, gff_example), store_sequences=True)

# scaffs = example_gff.scaffolds
# print(scaffs[0].seq[:50])
            
correct_n = 0
incorrect_n = 0
for seq_col in panaroo_obj.seq_collections:
    for gene in seq_col.sequences:
        gene_gff_file = os.path.join(gffs_dir, gene.get_genome_name()+".gff")
        gff_obj = parse_gff(gene_gff_file, store_sequences=True)
        # compare: cds seq, seq from scaffold splicing, seq from alignment
        # cds_sequences = gff_obj.get_sequences()[gene.seq_name]
        # cds_seq = gff_obj.get_cds_sequence(gene.annotation_id)
        scaff = gff_obj.get_scaffolds_dict()[gene.seq_name]
        csv_gene_seq = csv_dict.get(gene.annotation_id, "").seq
        csv_gene_coords = csv_dict.get(gene.annotation_id, "").refound_coords
        gene.seq = gene.seq.replace("-", "")
        # print("refound", gene.refound)
        # print("gff coords", gff_obj)
        # print("scaff_len", len(scaff.seq))
        # print("csv coords", csv_gene_coords)
        # # print("scaff", scaff.seq[:30])
        # print(gene.annotation_id, gene.start, gene.end, gene.strand)
        # print("panaroo_aln", gene.seq)
        # print("csv", csv_gene_seq)
        # # scaff_gene_seq = scaff.seq[gene.start:gene.end]
        
        if gene.strand < 0:
            scaff_gene_seq = Seq(scaff.seq).reverse_complement()
        else:
            scaff_gene_seq = Seq(scaff.seq)
        scaff_gene_seq = scaff_gene_seq[gene.start:gene.end]
        # print("scaff_spliced", scaff_gene_seq)
        if not (gene.seq.lower() == csv_gene_seq.lower() == scaff_gene_seq.lower()):
            incorrect_n += 1
            cds = gff_obj.get_cds_dict().get(gene.annotation_id, "NOT IN GFF")
            print(gene.annotation_id)
            print("refound", gene.refound)
            if not gene.refound:
                print("cds_coords", cds.start, cds.end)
                print("cds strand", cds.strand)
            print("gene_coords", gene.start, gene.end)
            print("strand", gene.strand)
            print("panaroo_aln\n", gene.seq.lower())
            print("csv\n", csv_gene_seq.lower())
            print("scaff\n", scaff_gene_seq.lower())
            print("="*30)
        else:
            correct_n += 1

print("corr, incorr", correct_n, incorrect_n)