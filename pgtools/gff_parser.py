import os
from Bio import SeqIO
from pgtools.utils import intersection_len, strand_rep
from typing import Dict

class GffCDS:
    def __init__(self, scaffold, annotation_id, strand, start, end):
        self.scaffold = scaffold
        self.ID = annotation_id
        self.strand = strand
        self.start = start
        self.end = end
    
    def __len__(self):
        return self.end - self.start + 1

class Scaffold:
    def __init__(self, genome, name, length):
        self.genome = genome
        self.name = name
        self.length = length
        self.seq = None
    
    def set_seq(self, seq: str):
        self.seq = seq

class SimpleAnnotation:
    def __init__(self, start, end, strand, annotation_id) -> None:
        self.start: int = start
        self.end: int = end
        self.strand: int = strand
        self.annotation_id: str = annotation_id
    
    def __str__(self):
        # strand_print = strand_rep(self.strand)
        # return strand_print
        str_strand = "+" if self.strand > 0 else "-"
        return f"{self.start} - {self.end}; strand: {str_strand}; {self.annotation_id}"
    
        

class Simple_GFF:
    def __init__(self, genome, scaffolds, CDSs):
        self.genome = genome
        scaffolds_coords = {scaff.name: set() for scaff in scaffolds}
        # assumes that in the same scaffold there are no cds with the same coords
        # as it is unlikely and allows for storing them in a set
        for cds in CDSs:
            scaffolds_coords[cds.scaffold.name].add(SimpleAnnotation(cds.start, cds.end, cds.strand, cds.ID))
        self.scaffolds_coords: Dict[str, SimpleAnnotation] = scaffolds_coords        

    def get_overlapping_annotations(self, return_pairs = False):
        gff = self.scaffolds_coords
        # scaff_cds = {scaff:{cds[-1]:[] for cds in coords} for scaff, coords in gff.items()}
        scaff_cds = {}
        for scaff, coords in gff.items():
            for cds_coords in coords:
                s = cds_coords[0]
                e = cds_coords[1]
                ann_id = cds_coords[-1]
                for cds in coords:
                    if intersection_len((s,e), (cds[0], cds[1])) > 0 and ann_id != cds[-1]:
                        if scaff not in scaff_cds:
                            scaff_cds[scaff] = {}
                        if ann_id not in scaff_cds[scaff]:
                            scaff_cds[scaff][ann_id] = []
                        scaff_cds[scaff][ann_id].append((cds[-1], intersection_len((s,e), (cds[0], cds[1]))))
        if return_pairs:
            #return dictionary saffold: {(cds1, cds2, overlap)}. No redundant information (but 
            # the default format comes in handy for comparisons)
            pass
        return scaff_cds        
    
    def get_highest_overlap(self):
        scaff_overlaps = self.get_overlapping_annotations()
        max_overlap = 0
        scaffold = None
        id_1 = None
        id_2 = None
        for scaff, cds in scaff_overlaps.items():
            for cds_, overlaps in cds.items():
                for el in overlaps:
                    if el[-1] > max_overlap:
                        max_overlap = el[-1]
                        id_1 = cds_
                        id_2 = el[0]
                        scaffold = scaff
        return {"scaffold":scaffold, "cds_1":id_1, "cds_2":id_2, "overlap":max_overlap}
    
class Pangenome_Gffs:
    def __init__(self, simple_GFFs):
        genome_gff = {}
        for gff in simple_GFFs:
            genome_gff[gff.genome] = gff
        self.simple_GFFs: Dict[str,Simple_GFF] = genome_gff
    
    def add_gff(self, Simple_GFF):
        self.simple_GFFs[Simple_GFF.genome] = Simple_GFF.scaffolds_coords
    
    def get_genome_cds_coords_dict(self, sep = ".") -> Dict[str, SimpleAnnotation]:
        all_cds_coords = {}
        for genome, gff in self.simple_GFFs.items():
            for scaff, cds in gff.scaffolds_coords.items():
                all_cds_coords[f"{genome}.{scaff}"] = set()
                for cds_ in cds:
                    all_cds_coords[f"{genome}.{scaff}"].add(cds_)
        return all_cds_coords

class Gff:
    """
    Represents gff file content
    """
    def __init__(self, scaffolds, cds):
        genomes = set([scaff.genome for scaff in scaffolds])
        assert len(genomes)==1, f"Provided scaffolds belong to more than one genome. Genomes: {genomes}"
        self.genome = genomes.pop()
        self.scaffolds = scaffolds
        self.cds = cds
    
    def get_scaffolds_dict(self):
        return {f"{scaffold.genome}.{scaffold.name}": scaffold for scaffold in self.scaffolds}
    
    def get_cds_dict(self):
        return {cds.ID: cds for cds in self.cds}

    def get_sequences(self):
        """
        Returns a dictionary Dict[genome.scaffold,cds]
        """
        sequences = {}

        for cds in self.cds:
            seq_name = f"{cds.scaffold.genome}.{cds.scaffold.name}"
            if not seq_name in sequences:
                sequences[seq_name] = [cds]
            else:
                sequences[seq_name].append(cds)
        
        return sequences

    def get_cds_sequence(self, gene_annotation_id: str) -> GffCDS:
        """
        NOT TESTED - are splicing coords coorect?
        ADD for - strand - reverse complement
        """
        for gene in self.cds:
            if gene.ID == gene_annotation_id:
                if gene.scaffold.seq:
                    seq = gene.scaffold.seq[gene.start + 1, gene.end]
                    return seq
        return None

    def join_gffs(self, other_gff):
        """
        Joins gff files
        """
        # ===TO DO======
        # check, if for the same genomes there are no conflicting cds and scaffold
        # ==============
    
    def to_MAF():
        """
        Writes cds into maf format
        Each seq is one block - does not make sense (no real MSA)
        but useful for filtering
        """
        

def strand_rep(strand_sign):
    if strand_sign == "+":
        return 1
    else:
        return -1

def parse_gff(gff_path, store_sequences=False):
    """
    Parses gff file

    Parameters:
        gff_path: str
            path to the gff file
    
    Returns:
        Gff
        Object containing information on scaffolds and annotations in gff file
    """
    genome_name = gff_path.strip("/")[-1].strip(".")[0]
    scaffolds_dict = {}
    CDS = []
    seq_block = False
    sequence = ""
    fasta_id = None
    with open(gff_path) as f:
        next(f)
        genome_name = gff_path.split("/")[-1][:-4]
        for line in f:

            if line.startswith("##FASTA"):
                if not store_sequences:
                    return Gff(list(scaffolds_dict.values()),CDS)            
                else:
                    seq_block = True
                    sequence = ""
            elif line.startswith("##sequence-region"):
                _, scaffold_name, _, scaffold_len = line.split()
                scaffolds_dict[scaffold_name] = Scaffold(genome_name, scaffold_name, int(scaffold_len))

            else:
                if not seq_block:
                    scaffold_name, _, _, start, end, _, strand, _, gene_info, *_ =line.split()
                    annotation_id = gene_info.split(";")[0][3:]

                    CDS.append(GffCDS(scaffolds_dict[scaffold_name], annotation_id, strand_rep(strand), int(start), int(end)))
                else:
                    if line.startswith(">"):
                        # print(line)
                        if sequence:
                            # print("seq",sequence[:50])
                            scaffolds_dict[fasta_id].seq = sequence
                            sequence = ""
                        fasta_id = line[1:].strip()
                    else:
                        sequence += line.strip()
    if sequence:
        scaffolds_dict[fasta_id].seq = sequence
    return Gff(list(scaffolds_dict.values()),CDS)            

def parse_GFFs_dir(GFFs_dir, gff_simple = True):
    """
    parses al gff files in a given directory
    """
    all_gff = {}
    for f in os.listdir(GFFs_dir):
        if f.endswith(".gff"):
            gff_obj = parse_gff(os.path.join(GFFs_dir, f), store_sequences=True)
            if gff_simple:
                all_gff[gff_obj.genome] = Simple_GFF(gff_obj.genome, gff_obj.scaffolds, gff_obj.cds)
            else:
                all_gff[gff_obj.genome] = gff_obj
    if gff_simple:
        return Pangenome_Gffs(all_gff.values())
    return all_gff


def scaffolds_from_GFFs_dir(GFFs_dir):
    """
    parses all gff files in a given directory
    """
    all_gff = {}
    for f in os.listdir(GFFs_dir):
        if f.endswith(".gff"):
            gff_obj = parse_gff(os.path.join(GFFs_dir, f), store_sequences=True)
            all_gff[gff_obj.genome] = gff_obj.scaffolds
    return all_gff