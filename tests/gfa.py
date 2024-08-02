from pgtools import gfa_parser
from pgtools import maf_parser
gfa_path = "/home/pampuch/studia/magisterka/pangenome_tools/data/tiny.gfa"
blocks, s_lines = gfa_parser.gfa1_to_maf(gfa_path)

# s_ids = [s.id for s in s_lines[1:]]
# s_unique_ids = list(set(s_ids))

# print(len(s_ids), len(s_unique_ids))

# print(len(blocks))

# gfa_v = []
# for i in range(1, len(blocks)+1):
#     gfa_v.append(gfa_parser.GfaVertex(s_lines[i], blocks[i]))

# gfa = gfa_parser.SimpleVertices(gfa_v)

gfa = gfa_parser.parse_gfa1(gfa_path)

genome_1 = blocks[1][0][0]
genome_2 = blocks[3][1][0]
print(genome_1, genome_2)

# g1_filtered = gfa.filter(genome_1)
# g2_filtered = gfa.filter(genome_2)

# print(g1_filtered)
# print(g2_filtered)
# common_v = g1_filtered.intersection(g2_filtered)
# print(common_v)

common_v = gfa.filter([genome_1, genome_2])

for id, V in common_v.items():
    V = gfa.vertices[id]
    # print(id)
    print(V.id)
    print("strand", V.seq_dict[genome_1].strand)

print("maf")
maf_path = "/home/pampuch/studia/magisterka/pangenome_tools/data/subset.maf"
maf = maf_parser.parse_maf(maf_path, store_seqs=False)
example_block = maf.synteny_blocks[3]
contig_names = [example_block.block_seqs[0].seq_name, example_block.block_seqs[0].seq_name]
maf_filtered = maf.filter(contig_names)

print(maf_filtered.synteny_blocks[0].id)
# maf_V = []
# for block in maf.synteny_blocks:
#     print(block)
#     print(block.id)
#     for seq in block.block_seqs:
#         s_line = gfa_parser.S_line(block.id, "", len(seq)
#         maf_V.append(gfa_parser.GfaVertex(gfa_parser.S_line(block)))

# maf_gfa = None