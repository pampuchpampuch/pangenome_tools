from itertools import product

class S_line:
	def __init__(self, name, sequence):
		self.id = name
		self.sequence = sequence
		self.length = len(sequence)

from collections import defaultdict

def gfa1_to_maf(gfa_file):
    s_lines = [None]
    blocks = defaultdict(list)
    srcSizes = defaultdict(int)
    with open(gfa_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("S"):
#				print "S line"
                line = line.split()
                name, sequence = line[1], line[2]
                s_lines.append(S_line(name, sequence))
            elif line.startswith("P"):
                line = line.split()
                alt = line[1]
                path = line[2].split(",")
                srcSize = 0
                for el in path:
                    vtx = int(el[:-1])
                    srcSize += s_lines[vtx].length
                srcSizes[alt] = srcSize
                pos = 0
                for el in path:
                    vtx = int(el[:-1])
                    strand = el[-1]
                    if strand == "-":
                        start = srcSize - pos - s_lines[vtx].length
                    else:
                        start = pos
                    pos += s_lines[vtx].length
                    blocks[vtx].append((alt, start, strand))
    return(blocks, s_lines)

def parse_gfa1(gfa_file):
    """
    From Ania's gfa1_to_maf script
    """
    s_lines = [None]
    blocks = defaultdict(list)
    srcSizes = defaultdict(int)
    with open(gfa_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("S"):
#				print "S line"
                line = line.split()
                name, sequence = line[1], line[2]
                s_lines.append(S_line(name, sequence))
            elif line.startswith("P"):
                line = line.split()
                alt = line[1]
                path = line[2].split(",")
                srcSize = 0
                for el in path:
                    vtx = int(el[:-1])
                    srcSize += s_lines[vtx].length
                srcSizes[alt] = srcSize
                pos = 0
                for el in path:
                    vtx = int(el[:-1])
                    strand = el[-1]
                    if strand == "-":
                        start = srcSize - pos - s_lines[vtx].length
                    else:
                        start = pos
                    pos += s_lines[vtx].length
                    blocks[vtx].append((alt, start, strand))
    # return(blocks, s_lines)
    gfa_v = []
    for i in range(1, len(blocks)+1):
        s_line, block = s_lines[i], blocks[i]
        id = int(s_line.id)
        length = s_line.length
        seq_dict = {}
        for seq in block:
            strand = 1 if seq[-1] == "+" else -1
            seq_dict[seq[0]] = SimpleSeq(seq[0], seq[1], seq[1] + length - 1, strand)
        gfa_v.append(GfaVertex(id, length, seq_dict))  

    return SimpleVertices(gfa_v)            

class SimpleSeq:
    def __init__(self, genome, start, end, strand):
        self.genome = genome
        self.start = start
        self.end = end
        self.strand = strand
       
class GfaVertex:
    # def __init__(self, s_line, block):
        # self.id = int(s_line.id)
        # self.length = s_line.length
        # seq_dict = {}
        # for seq in block:
        #     strand = 1 if seq[-1] == "+" else -1
        #     seq_dict[seq[0]] = SimpleSeq(seq[0], seq[1], seq[1] + self.length - 1, strand)
    def __init__(self, id, length, seq_dict):
        self.id = id
        self.length = length
        self.seq_dict = seq_dict

    def filter(self, seq_names):
        seq_dict = {seq_name: self.seq_dict[seq_name] for seq_name in seq_names
                    if seq_name in self.seq_dict.keys()}
        return GfaVertex(self.id, self.length, seq_dict)

    def __print__(self):
        print(f"{self.length, self.seq_dict}")

class SimpleVertices:
    def __init__(self, vertices):
        self.vertices = {V.id: V for V in vertices}

    def filter(self, genome_names):
        genome_names = set(genome_names)
        v_ids = set()
        filtered_v = []
        for id, V in self.vertices.items():
            # print(genome_names)
            # print(set(V.seq_dict.keys()))
            if genome_names.issubset(set(V.seq_dict.keys())):
                v_ids.add(id)
                filtered_v.append(V.filter(genome_names))
        return SimpleVertices(filtered_v)

    def get_filtered_vertices_by_strand(self,contig_names):
        gfa_verticles = {"".join(strand_comb): [] for strand_comb in product(["+","-"], repeat=len(contig_names))}
        for V in self.vertices.values():
            if not set(contig_names).issubset(set(V.seq_dict.keys())):
                continue
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
        
        return gfa_verticles