import os
import argparse
import networkx as nx
from Bio import SeqIO
import typing
from pgtools.maf_parser import MAFseq
from pgtools.gff_parser import GffCDS, Scaffold
from pgtools.pangenome import BaseSeq, SeqCollection, Pangenome

### TO DO ####
# 1. Object that keeps panaroo genes with information about start and end
# 2. Parse panaroo into a class, then tranform it into MAF object
#    - design class for panaroo output
#    - redesign class for PanarooCsvInfo
# 3. Parse more of panaroo output
##############
class PanarooCsvInfo(BaseSeq):
    """
    Represent gene information defined in gene_data.csv in panaroo output
    """
    def __init__(self, gff_file, scaffold_name, clustering_id, annotation_id, refound_coords=None, refound_strand=None, seq = None):
        self.gff_file = gff_file
        self.scaffold_name = scaffold_name
        self.chr_name = gff_file + "." + scaffold_name
        self.clustering_id = clustering_id
        self.annotation_id = annotation_id
        # if refound strand is none - gene is not refound (part of panaroo algorithm
        # detects ommited annotations)
        self.refound_strand = refound_strand
        self.refound_coords = refound_coords
        self.seq = seq

class PanarooGene(BaseSeq):
    """
    Represent gene as defined in panaroo model (gene_data.csv and other files in panaroo output)
    """

    def __init__(self, seq_name: str, start: int, end: int, strand: int, src_size: int, is_refound: bool, annotation_id: str, seq: str = None, in_format="panaroo", in_soft_core: bool = None):
        super().__init__(seq_name, start, end, strand, src_size, seq = seq, in_format = in_format, in_soft_core=in_soft_core)
        self.refound: bool = is_refound
        self.annotation_ids.append(annotation_id)
        self.seq = seq

class PanarooCluster(SeqCollection):
    """
    Represents panaroo cluster based on 
    """
    # def __init__(self, cluster_name, genes):
    #     self.name = cluster_name
    #     self.genes = genes

    # def __init__(self, id: int, cluster_name, seq_dict):
    #     # seq_dict = {seq.seq_name: seq for seq in genes}
    #     super().__init__(id, seq_dict)
    def __init__(self, id: int, seq_list: list[BaseSeq], cluster_name: str, is_soft_core: bool = None):
        super().__init__(id, seq_list, is_soft_core= is_soft_core)
        self.cluster_name: str = cluster_name
        # self.soft_core = is_soft_core
    def add_alignment(self, fasta):
        """
        Adds aligned sequences based on appropriate file in panaroo output
        aligned_gene_sequences folder
        """
        assert self.name == fasta.split("/")[0][:-4]

    def MAF_repr(self):
        pass

class Panaroo(Pangenome):
    """
    Represents Panaroo output
    """
    def __init__(self, seqs_collections: list[PanarooCluster], graph_structure: nx.Graph, soft_core_thresholds: typing.Dict[str, int]):
        super().__init__(seqs_collections)
        self.graph: nx.Graph = graph_structure
        self.soft_core_thresholds: typing.Dict[str, int] = soft_core_thresholds
    
    def get_clust_id_name_mapping(self):
        clust_id_mapping = {}
        for seq_coll in self.seq_collections:
            clust_id_mapping[seq_coll.id] = seq_coll.cluster_name

        return clust_id_mapping
    
    def assigned_annotations_to_csv(pg, gff_dir, csv_seqs_name = "panaroo_sequences_summary.csv", csv_annots_name = "panaroo_annotations_summary.csv"):
        res_csv = open(csv_seqs_name, "w")
        res_csv.write("cluster name,cluster size,core status,seq name,seq name,seq start,seq end,seq strand,CDS\n")
        res_csv.flush()
        annots_csv = open(csv_annots_name, "w")
        annots_csv.write("seq name,annotation ID\n")
        annots_csv.flush()
        pg.detect_soft_core()
        # pg.map_to_gff(gff_dir)
        for seq_coll in pg.seq_collections:
            id = seq_coll.cluster_name
            clust_size = len(seq_coll)
            core_status = seq_coll.soft_core
            for seq in seq_coll.sequences:
                # annotation len is also an iimportant aspect, but can be easily retrived from gff file
                # print(seq.seq_name, [ann.annotation_id for ann in seq.mapped_annotations], seq_coll.soft_core)
                annots = ";".join([ann_id for ann_id in seq.annotation_ids])
                res_csv.write(f"{id},{clust_size},{core_status},{seq.seq_name},{seq.start},{seq.end},{seq.strand},{annots}\n")
                res_csv.flush()
                for annot in seq.annotation_ids:
                    annots_csv.write(f"{seq.seq_name},{annot}\n")
                    annots_csv.flush()
        res_csv.close()
    
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

def parse_core_thresholds(panaroo_dir: str):
    core_thresholds = {}
    with open(os.path.join(panaroo_dir, "summary_statistics.txt")) as f:
        for line in f:
            category, thr_val = line.split("(")
            category = category.strip()
            thr_val = int(thr_val.split("%")[0])/100
            core_thresholds[category] = thr_val

    return core_thresholds

def parse_soft_core(panaroo_dir: str):
    file = "core_alignment_filtered_header.embl"
    core_clst_names = set()
    with open(os.path.join(panaroo_dir, file)) as f:
        clst_name = None
        for line in f:
            if not line.startswith("FT"):
                continue
            if line.split()[1] == "feature":
                if clst_name:
                    core_clst_names.add(clst_name)
            else:
                if line.split()[1].startswith("/label"):
                    clst_name = line.split("=")[-1].strip()[:-len(".aln")]
        core_clst_names.add(clst_name)
            
    return core_clst_names

def parse_gene_data(gene_data_path):
    """
    Parses gene_data.csv (part of Panaroo output)
    """
    genes_dict = {}

    with open(gene_data_path) as f:
        next(f)
        for line in f:
            if line:
                # gff_file, scaffold_name, clustering_id, annotation_id, _, _, _, description, additional_info = line.split(",")
                line_items = line.split(",")
                gff_file, scaffold_name, clustering_id, annotation_id, _, seq, *_ = line.split(",")
                
                refound_coords = None
                refound_strand = None
                if clustering_id.split("_")[1] == "refound":
                    description, additional_info = line_items[-1].split(";")
                    refound_coords = tuple(map(lambda x: int(x), description[len("location:"):].split("-")))
                    refound_strand = strand_rep(additional_info[len("strand:")])
                genes_dict[clustering_id] = PanarooCsvInfo(gff_file, scaffold_name, clustering_id, annotation_id, refound_coords, refound_strand, seq=seq)

    return genes_dict

def parse_panaroo_aln(aln_file: str, gff_dict: typing.Dict[str, GffCDS], scaffolds_dict: typing.Dict[str, Scaffold], genes_dict: typing.Dict[str, PanarooCsvInfo], cluster_id: int, include_refound = True, set_soft_core:typing.Set[str] = False) -> PanarooCluster:
    """
    Parses one Panaroo cluster (based on output files
    in the panaroo aligned_gene_sequences)
    """
    panaroo_genes = []
    seq_lens = []

    clst_name  = aln_file.split("/")[-1].split(".")[0]
    is_soft_core = None
    if set_soft_core:
        is_soft_core = True if clst_name in set_soft_core else False
    # print(clst_name)
    # print(set_soft_core)
    # print(is_soft_core)
    with open(aln_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_name, clustering_id = record.id.split(";")
            seq = record.seq
            seq_lens.append(len(seq))
            if seq_name[:3] == "_R_":
                # print("R")
                strand = -1
                seq_name = seq_name[3:]
            else:
                strand = 1

            gene_info  = genes_dict[clustering_id]
            # try:
            #     cds_info = gff_dict[gene_info.annotation_id]
            # except:
            #     # print(f"Sequence with annotation id {gene_info.annotation_id} from alignment file {aln_file} not found in gff files, ommiting the sequence from block")
            #     continue

            is_refound = False
            if gene_info.refound_strand:
                if include_refound:
                ###  HOW TO GET CONTIG SIZE
                    is_refound = True
                    chr_name = gene_info.gff_file + "." + gene_info.scaffold_name
                    chr_size = scaffolds_dict[chr_name].length
                    start, end = gene_info.refound_coords
                    # refound coord system is not in agreement with other coords in panaroo
                    # is zero based and if strand is "-", coords are given for rev strand
                    strand = gene_info.refound_strand * strand
                # if store_annotation:
                #     panaroo_genes.append(PanarooGene(chr_name + ";" + gene_info.annotation_id, start, end, strand, chr_size, in_format="panaroo_refound", is_refound=is_refound, annotation_id=gene_info.annotation_id, seq = seq, in_soft_core=is_soft_core))

                # else:
                panaroo_genes.append(PanarooGene(chr_name, start, end, strand, chr_size, in_format="panaroo_refound", is_refound=is_refound, annotation_id=gene_info.annotation_id, seq = seq, in_soft_core=is_soft_core))
            else:
                cds_info = gff_dict[gene_info.annotation_id]   
                chr_name = gene_info.chr_name
                chr_size = cds_info.scaffold.length
                strand = cds_info.strand * strand
                start = cds_info.start
                end = cds_info.end
                # if store_annotation:
                #     panaroo_genes.append(PanarooGene(chr_name + ";" + gene_info.annotation_id, start, end, strand, chr_size, in_format="panaroo", is_refound=is_refound, annotation_id=gene_info.annotation_id, seq = seq, in_soft_core=is_soft_core))

                # else:
                panaroo_genes.append(PanarooGene(chr_name, start, end, strand, chr_size, in_format="panaroo", is_refound=is_refound, annotation_id=gene_info.annotation_id, seq = seq, in_soft_core=is_soft_core))
            # if strand > 0:
            #     start = start - 1
            #     end = end          
            # else:
            #     start = chr_size - end
            #     end = chr_size - start + 1

            # if cds_info.start != start:
                # print(chr_size)
            # print(f"cds: {cds_info.start}-{cds_info.end}")
            # print(f"maf: {start}-{end}")
            
            # seq_no_gaps = "".join(panaroo_gene.seq.split("-"))
            # msg = f"{panaroo_gene.start}:{panaroo_gene.end},{cds_info.start}:{cds_info.end}, {len(cds_info)}, {len(panaroo_gene)}, {len(seq_no_gaps)}"
            # assert len(panaroo_gene)==len(cds_info)==len(seq_no_gaps), print(msg)

            # # print(panaroo_gene.start,panaroo_gene.end)
            # # print(panaroo_gene.gff_coords(panaroo_gene.start, panaroo_gene.end, panaroo_gene.chr_size, panaroo_gene.strand))
            # # print(cds_info.start, cds_info.end)
            # maf_s, maf_e = panaroo_gene.gff_coords(panaroo_gene.start, panaroo_gene.end, panaroo_gene.chr_size, panaroo_gene.strand)
            # assert maf_s == cds_info.start and maf_e == cds_info.end, print(maf_s, cds_info.start)
                
            # panaroo_genes.append(panaroo_gene)
        # maf_out.write('a\n')
        # for maf_seq in maf_seqs:
        #     maf_out.write(maf_seq.MAF_repr())
        # maf_out.write('\n')
        first_seq_len = seq_lens[0]
        for seq_len in seq_lens:
            assert seq_len == first_seq_len, seq_lens
        # print("seq lens ok")

        return PanarooCluster(cluster_id, panaroo_genes, cluster_name = clst_name, is_soft_core = is_soft_core)

def parse_panaroo_output(panaroo_dir: str, gffs_dir: str, include_refound = True, detect_core = True) -> Panaroo:
    """
    Parses panaroo output
    """
    core_thresholds = parse_core_thresholds(panaroo_dir)
    # print(core_thresholds)
    soft_core_cluster_names = parse_soft_core(panaroo_dir)
    # print(soft_core_cluster_names)
    gff_dict = {}
    scaffolds_dict = {}
    for filename in os.listdir(gffs_dir):
        scaff_dict_, gff_dict_ = parse_gff(os.path.join(gffs_dir,filename))
        gff_dict.update(gff_dict_)
        scaffolds_dict.update(scaff_dict_)

    gene_data_dir = os.path.join(panaroo_dir, "gene_data.csv")

    genes_dict = parse_gene_data(gene_data_dir)

    aln_files_dir = os.path.join(panaroo_dir, "aligned_gene_sequences")

    panaroo_clusters = []

    clust_id = 0
    for aln_file in os.listdir(aln_files_dir):
        aln_file = os.path.join(aln_files_dir,aln_file)
        if not detect_core:
            pan_clust = parse_panaroo_aln(aln_file, gff_dict, scaffolds_dict, genes_dict, clust_id, include_refound=include_refound)
        else:
            pan_clust = parse_panaroo_aln(aln_file, gff_dict, scaffolds_dict, genes_dict, clust_id, include_refound=include_refound, set_soft_core=soft_core_cluster_names)
        panaroo_clusters.append(pan_clust)
        clust_id += 1
    
    G = nx.read_gml(os.path.join(panaroo_dir, "final_graph.gml"))
    return Panaroo(panaroo_clusters, G, core_thresholds)

def panaroo_aln_to_maf(aln_file: str, gff_dict: typing.Dict[str, GffCDS], scaffolds_dict: typing.Dict[str, Scaffold], genes_dict: typing.Dict[str, PanarooCsvInfo], maf_out: str, include_refound = True):
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
            # try:
            #     cds_info = gff_dict[gene_info.annotation_id]
            # except:
            #     # print(f"Sequence with annotation id {gene_info.annotation_id} from alignment file {aln_file} not found in gff files, ommiting the sequence from block")
            #     continue
            if gene_info.refound_strand and include_refound:
                ###  HOW TO GET CONTIG SIZE
                chr_name = gene_info.gff_file + "." + gene_info.scaffold_name
                chr_size = scaffolds_dict[chr_name].length
                start, end = gene_info.refound_coords
                strand = gene_info.refound_strand

                if strand > 0:
                    start = cds_info.start - 1
                    end = cds_info.end          
                else:
                    start = chr_size - cds_info.end
                    end = chr_size - cds_info.start + 1

                maf_seq = MAFseq(chr_name, start, end, strand, chr_size, str(seq.upper()))

            else:
                cds_info = gff_dict[gene_info.annotation_id]   
                chr_name = gene_info.chr_name
                chr_size = cds_info.scaffold.length

                # print(cds_info.strand,strand)

                strand = cds_info.strand * strand
                if strand > 0:
                    start = cds_info.start - 1
                    end = cds_info.end          
                else:
                    start = chr_size - cds_info.end
                    end = chr_size - cds_info.start + 1

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

def panaroo_output_to_maf(panaroo_dir, gff_dir, maf_out):
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
            


