import os
from Bio import SeqIO

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
