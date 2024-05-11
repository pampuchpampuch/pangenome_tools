import os
from Bio import SeqIO
import warnings

class Gff:
    def __init__(self, scaffolds, cds):
        genomes = set([scaff.genome for scaff in scaffolds])
        assert len(genomes)==1, f"Provided scaffolds belong to more than one genome. Genomes: {genomes}"
        self.genome = genomes.pop()
        self.scaffolds = scaffolds
        self.cds = cds
    
    def get_sequences(self):
        """
        Returns a dictioanry Dict[genome.scaffold,cds]
        """
        sequences = {}

        for cds in self.cds:
            seq_name = f"{cds.scaffold.genome}.{cds.scaffold.name}"
            if not seq_name in sequences:
                sequences[seq_name] = [cds]
            else:
                sequences[seq_name].append(cds)
        
        return sequences

    def join_gffs(self, other_gff):
        """
        Joins gff files
        """
        # ===TO DO======
        # check, if for the same genomes there are no conflicting cds and scaffold
        # ==============


class Scaffold:
    def __init__(self, genome, name, length):
        self.genome = genome
        self.name = name
        self.length = length

class GffCDS:
    def __init__(self, scaffold, annotation_id, strand, start, end):
        self.scaffold = scaffold
        self.ID = annotation_id
        self.strand = strand
        self.start = int(start)
        self.end = int(end)

def strand_rep(strand_sign):
    if strand_sign == "+":
        return 1
    else:
        return -1

def parse_gff(gff_path):
    scaffolds_dict = {}
    CDS = []

    with open(gff_path) as f:
        next(f)
        genome_name = gff_path.split("/")[-1][:-4]
        for line in f:

            if line.startswith("##FASTA"):
                return Gff(list(scaffolds_dict.values()),CDS)            
            if line.startswith("##sequence-region"):
                _, scaffold_name, _, scaffold_len = line.split()
                scaffolds_dict[scaffold_name] = Scaffold(genome_name, scaffold_name, scaffold_len)

            else:
                scaffold_name, _, _, start, end, _, strand, _, gene_info, *_ =line.split()
                annotation_id = gene_info.split(";")[0][3:]

                CDS.append(GffCDS(scaffolds_dict[scaffold_name], annotation_id, strand_rep(strand), start, end))

    return Gff(list(scaffolds_dict.values()),CDS)            
