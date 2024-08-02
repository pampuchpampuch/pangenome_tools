import os
import argparse
from Bio import SeqIO
from pgtools.maf_parser import MAFseq
from pgtools.gff_parser import GffCDS, Scaffold

### TO DO ####
# 1. Object that keeps panaroo genes with information about start and end
# 2. Parse panaroo into a class, then tranform it into MAF object
#    - design class for panaroo output
#    - redesign class for PanarooGene
# 3. Parse more of panaroo output
##############

class Panaroo:
    """
    Represents Panaroo output
    """
    pass

class PanarooGene:
    """
    Represent gene as defined in panaroo model (gene_data.csv in panaroo output)
    """
    def __init__(self, gff_file, scaffold_name, clustering_id, annotation_id, seq = None):
        self.gff_file = gff_file
        self.scaffold_name = scaffold_name
        self.chr_name = gff_file + "." + scaffold_name
        self.clustering_id = clustering_id
        self.annotation_id = annotation_id
        self.seq = seq

class PanarooCluster:
    """
    Represents panaroo cluster based on 
    """
    def __init__(self, cluster_name, genesm ):
        self.name = cluster_name
        self.genes = genes

    def add_alignment(self, fasta):
        """
        Adds aligned sequences based on appropriate file in panaroo output
        aligned_gene_sequences folder
        """
        assert self.name == fasta.split("/")[0][:-4]

    def MAF_repr(self):
        pass

def strand_rep(strand_sign):
    if strand_sign == "+":
        return 1
    else:
        return -1
    
def parse_gff(gff_path):
    """
    Parses gff file in a manner suitable for retrieving information
    from Panaroo output
    """
    scaffolds_dict = {}
    CDS_dict = {}

    with open(gff_path) as f:
        next(f)
        genome_name = gff_path.split("/")[-1][:-4]
        for line in f:

            if line.startswith("##FASTA"):
                return scaffolds_dict, CDS_dict
            
            if line.startswith("##sequence-region"):
                _, scaffold_name, _, scaffold_len = line.split()
                scaffolds_dict[genome_name+"."+scaffold_name] = Scaffold(genome_name, scaffold_name, int(scaffold_len))
            else:
                scaffold_name, _, _, start, end, _, strand, _, gene_info, *_ =line.split()
                annotation_id = gene_info.split(";")[0][3:]

                CDS_dict[annotation_id] = GffCDS(scaffolds_dict[genome_name+"."+scaffold_name], annotation_id, strand_rep(strand), int(start), int(end))

    return scaffolds_dict, CDS_dict

def parse_gene_data(gene_data_path):
    """
    Parses gene_data.csv (part of Panaroo output)
    """
    genes_dict = {}

    with open(gene_data_path) as f:
        next(f)
        for line in f:
            if line:
                gff_file, scaffold_name, clustering_id, annotation_id, *_ = line.split(",")
                genes_dict[clustering_id] = PanarooGene(gff_file, scaffold_name, clustering_id, annotation_id)
    
    return genes_dict

def panaroo_aln_to_maf(aln_file, gff_dict, genes_dict, maf_out):
    """
    Writes maf block based on one Panaroo cluster (based on output files
    in the panaroo aligned_gene_sequences)
    """
    maf_seqs = []
    with open(aln_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_name, clustering_id = record.id.split(";")
            seq = record.seq

            if seq_name[:3] == "_R_":
                # print("R")
                strand = -1
                seq_name = seq_name[3:]
            else:
                strand = 1

            gene_info  = genes_dict[clustering_id]
            try:
                cds_info = gff_dict[gene_info.annotation_id]
            except:
                # print(f"Sequence with annotation id {gene_info.annotation_id} from alignment file {aln_file} not found in gff files, ommiting the sequence from block")
                continue
            chr_name = gene_info.chr_name
            chr_size = cds_info.scaffold.length

            # print(cds_info.strand,strand)

            strand = cds_info.strand * strand
            if strand > 0:
                start = cds_info.start
                end = cds_info.end          
            else:
                start = chr_size - cds_info.end
                end = chr_size - cds_info.start

            # if cds_info.start != start:
                # print(chr_size)
            # print(f"cds: {cds_info.start}-{cds_info.end}")
            # print(f"maf: {start}-{end}")
            maf_seq = MAFseq(chr_name, start, end, strand, chr_size, str(seq.upper()))
            seq_no_gaps = "".join(maf_seq.seq.split("-"))
            msg = f"{maf_seq.start}:{maf_seq.end},{cds_info.start}:{cds_info.end}, {len(cds_info)}, {len(maf_seq)}, {len(seq_no_gaps)}"
            assert len(maf_seq)==len(cds_info)==len(seq_no_gaps), print(msg)

            # print(maf_seq.start,maf_seq.end)
            # print(maf_seq.gff_coords(maf_seq.start, maf_seq.end, maf_seq.chr_size, maf_seq.strand))
            # print(cds_info.start, cds_info.end)
            maf_s, maf_e = maf_seq.gff_coords(maf_seq.start, maf_seq.end, maf_seq.chr_size, maf_seq.strand)
            assert maf_s == cds_info.start and maf_e == cds_info.end, print(maf_s, cds_info.start)
            
            maf_seqs.append(maf_seq)
        
        maf_out.write('a\n')
        for maf_seq in maf_seqs:
            maf_out.write(maf_seq.MAF_repr())
        maf_out.write('\n')

def parse_panaroo_output(panaroo_dir, gff_dir, maf_out):
    """
    Writes maf file based on panaroo output
    """
    gff_dict = {}
    for filename in os.listdir(gff_dir): gff_dict.update(parse_gff(os.path.join(gff_dir,filename))[1])
    gene_data_dir = os.path.join(panaroo_dir, "gene_data.csv")
    genes_dict = parse_gene_data(gene_data_dir)
    aln_files_dir = os.path.join(panaroo_dir, "aligned_gene_sequences")

    maf_out = open(maf_out,"w")
    maf_out.write("##maf version=1 scoring=N/A\n\n")

    for aln_file in os.listdir(aln_files_dir):
        aln_file = os.path.join(aln_files_dir,aln_file)
        panaroo_aln_to_maf(aln_file, gff_dict, genes_dict, maf_out)

    maf_out.close()

def get_vertices():
    pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("panaroo_dir",
                        help="path to dir with panaroo output")
    parser.add_argument("gff_dir",
                        help="path to dir with gff files used to generate panaroo output")
    parser.add_argument("maf_out",
                    help="path to output maf file")
    args = parser.parse_args()

    gff_dict = {}
    for filename in os.listdir(args.gff_dir): gff_dict.update(parse_gff(os.path.join(args.gff_dir,filename))[1])

    gene_data_dir = os.path.join(args.panaroo_dir, "gene_data.csv")

    genes_dict = parse_gene_data(gene_data_dir)

    aln_files_dir = os.path.join(args.panaroo_dir, "aligned_gene_sequences")

    maf_out = open(args.maf_out,"w")
    maf_out.write("##maf version=1 scoring=N/A\n\n")

    for aln_file in os.listdir(aln_files_dir):
        aln_file = os.path.join(aln_files_dir,aln_file)
        panaroo_aln_to_maf(aln_file, gff_dict, genes_dict, maf_out)

    maf_out.close()


if __name__ == "__main__":
    main()
            


