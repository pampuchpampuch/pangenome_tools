from pgtools import gfa_parser
from pgtools import maf_parser

def intersection_len(s1_coords, s2_coords):
    """
    seq2 is assumed to be from MAF and sed1 from gfa (checking how well
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

def contains(s1_coords, s2_coords):
    """
    Assumes s1 is from gfa and s2 from maf. The goal is to check if
    gfa vertex is contained in maf block
    """
    s1_start, s1_end = s1_coords
    s2_start, s2_end = s2_coords
    return (s1_start >= s2_start) and (s1_end <= s2_end)

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

gfa_filtered = gfa.filter([genome_1, genome_2])

for id, V in gfa_filtered.vertices.items():
    V = gfa.vertices[id]
    # print(id)
    print(V.id)
    print("strand", V.seq_dict[genome_1].strand)

print("maf")
maf_path = "/home/pampuch/studia/magisterka/pangenome_tools/data/subset.maf"
maf = maf_parser.parse_maf(maf_path, store_seqs=False)

example_block = maf.synteny_blocks[0]
example_block_size  = len(example_block.block_seqs)

for block in maf.synteny_blocks:
    if len(block.block_seqs) > example_block_size:
        example_block = block
        example_block_size = len(block.block_seqs)

n_seqs = len(example_block.block_seqs)
all_contig_names = [example_block.block_seqs[i].seq_name for i in range(n_seqs)]
print(all_contig_names)

#change names in gfa verticles and maf blocks so that it is easier to find matching ones

verticles_alt = []


for id, V in gfa.vertices.items():
    alt_seq_dict = {}
    if not id % 3 == 0:
        continue
    i = 0
    for genome_name, seq in V.seq_dict.items():
        alt_seq = seq
        alt_seq.genome = all_contig_names[i]
        alt_seq_dict[all_contig_names[i]] = alt_seq
        i += 1
    verticles_alt.append(gfa_parser.GfaVertex(id, V.length, alt_seq_dict))

# print(verticles_alt)
for V in verticles_alt:
    for genome, seq in V.seq_dict.items():
        print(genome)
        print(seq.genome)

gfa_alt = gfa_parser.SimpleVertices(verticles_alt)
contig_names = all_contig_names[:2]

gfa_alt_filtered = gfa_alt.filter(contig_names)
maf_filtered = maf.filter(contig_names)

print(f"gfa before: {len(gfa.vertices)}, after: {len(gfa_alt_filtered.vertices)}")
for block in maf.synteny_blocks[1:3]:
    for i in range(len(block.block_seqs)):
        block.block_seqs[i].seq_name = all_contig_names[i]

# contig_names = all_contig_names[:2]

print(maf_filtered.synteny_blocks[0].id)
# maf_V = []
for block in maf_filtered.synteny_blocks:
    print(block)
    print(block.id)
    # for seq in block.block_seqs:
        # s_line = gfa_parser.S_line(block.id, "", len(seq)
        # maf_V.append(gfa_parser.GfaVertex(gfa_parser.S_line(block)))

# maf_gfa = None
print(f"n blocks before {len(maf.synteny_blocks)}, after: {len(maf_filtered.synteny_blocks)}")

# would be best to check sorted: if one contig does not fit in first k blocks, the next one
# will not fit either

"""
blocks: categorise by strands:
    blocks ++. list of pairs of coords [((s1, e1),(s2, e2))]
    blocks +-
    blocks -+
    blocks --

The same for gfa

"""
gfa_verticles = {"++": [], "+-": [], "--": [], "-+": []}
for V in gfa_alt_filtered.vertices.values():
    dict_key = ""
    coords = []
    for genome in contig_names:
        seq = V.seq_dict[genome]
        coords.append((seq.start, seq.end))
        if seq.strand > 0:
            dict_key += "+"
        else:
            dict_key += "-"

    gfa_verticles[dict_key].append(tuple(coords))

print(gfa_verticles)

maf_blocks = {"++": [], "+-": [], "--": [], "-+": []}

for block in maf_filtered.synteny_blocks:
    dict_key = ""
    coords = []
    seqs_dict = {seq.seq_name: seq for seq in block.block_seqs}
    # print(seqs_dict)
    for genome in contig_names:
        seq = seqs_dict[genome]
        coords.append((seq.start, seq.end))
        if seq.strand > 0:
            dict_key += "+"
        else:
            dict_key += "-"       

    maf_blocks[dict_key].append(tuple(coords))
print("MAF")
print(maf_blocks)

maf_mock = [((100,500),(110,550)), ((1000,2000),(1010,2090))]
maf_mock_coords =  {"++": maf_mock,
                    "+-": [],
                    "--": [],
                    "-+": []}

gfa_mock_coords =  {"++": [((34,67),(145,490)), ((110,300),(120,400))],
                    "+-": [],
                    "--": [],
                    "-+": []}


contained_lens = []
not_contained_lens = []

def vertices_freq(maf_blocks, gfa_verticles):
    for strandness in maf_blocks.keys():
        for V in gfa_verticles[strandness]:
            in_block = False
            for block in maf_blocks[strandness]:
                print("v",V)
                print("block", block)
                if contains(V[0], block[0]) and contains(V[1], block[1]):
                    contained_lens.append(V[0][1]-V[0][0])
                    in_block = True
            if not in_block:
                not_contained_lens.append(V[0][1]-V[0][0])
    return contained_lens, not_contained_lens


cont, not_cont = vertices_freq(maf_mock_coords, gfa_mock_coords)
print(cont)
print(not_cont)
