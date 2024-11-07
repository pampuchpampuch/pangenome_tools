import typing
from itertools import product
import numpy as np
import copy
from pgtools.gff_parser import parse_GFFs_dir
from pgtools.utils import contains
import networkx as nx
from Bio import SeqIO
import os
import subprocess
import re

### OLD WORKING VERSION #####

FORMATS_COORD_SYSTEMS = {
    "maf": {"coord_system": ("0-based", "half-open"), "rev_strand_coords": "rev_strand"},
    # per ensembl, gff is NOT half open (but fully-closed) - but works???
    "gff": {"coord_system": ("1-based", "fully-closed"), "rev_strand_coords": "forward_strand"},
    "gfa": {"coord_system": ("0-based", "half-open"), "rev_strand_coords": "forward_strand"},
    # is panaroo coord system the same for refound and regular (gff) genes - NO, refound is 0-based?
    "panaroo": {"coord_system": ("1-based", "half-open"), "rev_strand_coords": "forward_strand"},
    "panaroo_refound":  {"coord_system": ("0-based", "half-open"), "rev_strand_coords": "forward_strand"}
    }

def invert_coords(start, end, seq_size, format):
    assert format in FORMATS_COORD_SYSTEMS, format
    coord_system = FORMATS_COORD_SYSTEMS[format]["coord_system"]
    # what if end is not half open?????
    # in theory, gff is fully closed, but Panaroo seems to
    # assume that it is not (or mabey it is because of how splicing works?)
    if coord_system[0] == "0-based":
        res_end = seq_size - start
    else:
        res_end = seq_size - start + 1
    if coord_system[1] == "fully-closed":
        res_start = seq_size - end + 1
    else:
        res_start = seq_size - end
    return res_start, res_end


def convert_coords(start: int, end: int, strand: int, src_size: int, in_format: str, out_format: str = "maf"):
    ### !!! other end than half-open is NOT supported!
    in_coord_system = FORMATS_COORD_SYSTEMS[in_format]["coord_system"]
    out_coord_system = FORMATS_COORD_SYSTEMS[out_format]["coord_system"]
    in_rev = FORMATS_COORD_SYSTEMS[in_format]["rev_strand_coords"]
    out_rev = FORMATS_COORD_SYSTEMS[out_format]["rev_strand_coords"]
    # if in is rev strand and the out is not: reverse if strand is -

    # start working on forward strand
    if in_rev == "rev_strand" and strand < 0:
        start, end = invert_coords(start, end, src_size, in_format)
              
    # then to "universal" coord system (0-based, half open, forward strand)

    if in_coord_system[0] == "1-based":
        start -= 1

    # then to chosen out format 

    if out_coord_system[0] == "1-based":
        start += 1
    if out_rev == "rev_strand" and strand < 0:
        start, end = invert_coords(start, end, src_size, out_format)

    return start, end

#### END OLD WORKING VERSION #####
def parse_mafft_output(alignment_fasta):
    '''
    parses mafft output - one MSA
    '''

    # seq_list=re.findall(">[^>]+",alignment)
    # sequences=[]

    # for seq in seq_list:
    #     id_end=re.match(".+",seq).end()
    #     sequence=seq[id_end:].replace('\n','')
    #     if sequence:
    #         sequences.append(sequence)

    # return(sequences)

    records = list(SeqIO.parse(alignment_fasta, "fasta"))
    return {record.id: record.seq for record in records}

class BaseSeq:
    """
    Represent a genomic sequence as a fragment of bigger Source fragment
    (usually chromosome or a contig)
    As the MAF files were commonly used throughout the analysis, default coordination system
    for this class is the one used in MAF files (zero-based, half open and coords for "-" strand
    are given in respect to the reversed strand)
    """
    def __init__(self, seq_name: str, start: int , end: int, strand: int, src_size: int, in_format, seq: str = None, coord_system: str ="maf", in_soft_core: bool = None):
        assert coord_system in FORMATS_COORD_SYSTEMS
        self.coord_system = coord_system
        self.seq_name: str = seq_name
        assert strand == -1 or strand == 1
        self.strand: int = strand
        start, end = convert_coords(start, end, self.strand, src_size, in_format=in_format)
        self.start: int = start
        self.end: int = end
        self.src_size: int = src_size
        self.seq: str = seq
        self.in_soft_core: bool = in_soft_core
        # description can store additional info ie mapped cds
        self.annotation_ids: typing.List[str] = []
        self.mapped_annotations: typing.List[str] = []
        self.cluster_id = None
        ### Add coord system?
        ### For which strands are coords given in the case of - strand

    def invert_coords(self):
        start = self.start
        end = self.end
        seq_size = self.src_size
        format = self.coord_system
        assert format in FORMATS_COORD_SYSTEMS, format
        coord_system = FORMATS_COORD_SYSTEMS[format]["coord_system"]
        # what if end is not half open?????
        # in theory, gff is fully closed, but Panaroo seems to
        # assume that it is not (or mabey it is because of how splicing works?)
        if coord_system[0] == "0-based":
            res_end = seq_size - start
        else:
            res_end = seq_size - start + 1
        if coord_system[1] == "fully-closed":
            res_start = seq_size - end + 1
        else:
            res_start = seq_size - end
        return res_start, res_end


    # # def convert_coords(start: int, end: int, strand: int, src_size: int, in_format: str, out_format: str = "maf"):
    def convert_coords(self, out_format):
        in_format = self.coord_system
        strand = self.strand
        start = self.start
        end = self.end
        src_size = self.src_size

        res_seq = copy.deepcopy(self)
        ### !!! other end than half-open is NOT supported!
        in_coord_system = FORMATS_COORD_SYSTEMS[in_format]["coord_system"]
        out_coord_system = FORMATS_COORD_SYSTEMS[out_format]["coord_system"]
        in_rev = FORMATS_COORD_SYSTEMS[in_format]["rev_strand_coords"]
        out_rev = FORMATS_COORD_SYSTEMS[out_format]["rev_strand_coords"]
        # if in is rev strand and the out is not: reverse if strand is -

        # start working on forward strand
        if in_rev == "rev_strand" and strand < 0:
            start, end = invert_coords(start, end, src_size, in_format)
                
        # then to "universal" coord system (0-based, half open, forward strand)

        if in_coord_system[0] == "1-based":
            start -= 1

        # then to chosen out format 

        if out_coord_system[0] == "1-based":
            start += 1
        if out_rev == "rev_strand" and strand < 0:
            start, end = invert_coords(start, end, src_size, out_format)

        res_seq.start = start
        res_seq.end = end
        res_seq.coord_system = out_format

        return res_seq
    # def invert_coords(self):
    #     assert format in FORMATS_COORD_SYSTEMS
    #     coord_system = FORMATS_COORD_SYSTEMS[format]["coord_system"]
    #     start = self.start
    #     end = self.end
    #     seq_size = self.src_size
    #     # what if end is not half open?????
    #     # in theory, gff is fully closed, but Panaroo seems to
    #     # assume that it is not (or mabey it is because of how splicing works?)
    #     if coord_system[0] == "0-based":
    #         res_end = seq_size - start
    #     else:
    #         res_end = seq_size - start + 1
    #     if coord_system[1] == "fully-closed":
    #         res_start = seq_size - end + 1
    #     else:
    #         res_start = seq_size - end
    #     return res_start, res_end

    # def convert_coords(self, out_format):
    #     in_format = self.coord_system
    #     ### !!! other end than half-open is NOT supported!
    #     in_coord_system = FORMATS_COORD_SYSTEMS[in_format]["coord_system"]
    #     out_coord_system = FORMATS_COORD_SYSTEMS[out_format]["coord_system"]
    #     in_rev = FORMATS_COORD_SYSTEMS[in_format]["rev_strand_coords"]
    #     out_rev = FORMATS_COORD_SYSTEMS[out_format]["rev_strand_coords"]
    #     # if in is rev strand and the out is not: reverse if strand is -
    #     strand = self.strand
    #     start = self.start
    #     end = self.end
    #     src_size = self.src_size
    #     # start working on forward strand
    #     if in_rev == "rev_strand" and strand < 0:
    #         start, end = invert_coords(start, end, src_size, in_format)
                
    #     # then to "universal" coord system (0-based, half open, forward strand)

    #     if in_coord_system[0] == "1-based":
    #         start -= 1

    #     if in_coord_system[1] == "fully-closed":
    #         end += 1

    #     # then to chosen out format 

    #     if out_coord_system[0] == "1-based":
    #         start += 1
    #     if out_coord_system[1] == "fully-closed":
    #         end -= 1
    #     if out_rev == "rev_strand" and strand < 0:
    #         start, end = invert_coords(start, end, src_size, out_format)

        return start, end    
    def get_genome_name(self, sep="."):
        return self.seq_name.split(sep)[0]
    
    def get_contig_name(self, sep = "."):
        return self.seq_name.split(sep)[1]

    def set_soft_core(self, is_soft_core: bool):
        self.in_soft_core = is_soft_core      
    
    def __str__(self):
        return (f"{self.seq_name}: {self.start} - {self.end}, {self.strand}")

    def __len__(self):
        """
        only half-open coord systems supported. For fully closed minor changes
        would have to made
        """
        # assert self.coord_system == "maf"
        format_coords = FORMATS_COORD_SYSTEMS[self.coord_system]["coord_system"]
        if format_coords[0] == "0-based":
            return self.end - self.start
        else:
            return self.end - self.start + 1

    def MAF_repr(self):
        """
        returns string representing seguence in MAF format
        """
        
        strand_sign = "+" if self.strand > 0 else "-"

        s_line=f"s\t{self.seq_name}\t{self.start}\t{len(self)}\t{strand_sign}\t{self.src_size}\t{self.seq}\n"

        return(s_line)
            

class SeqCollection:
    def __init__(self, id: int, seq_list: list[BaseSeq], is_soft_core: bool = None):
        self.id: typing.Any = id
        # seq_dict = {seq.seq_name: seq for seq in seq_list}
        seq_dict = {seq for seq in seq_list}        
        self.sequences: typing.Set[BaseSeq] = set(seq_list)
        # self.genome_names: typing.Set[str] = set(seq.get_genome_name() for seq in seq_dict.values())
        self.soft_core: bool = is_soft_core

    def __len__(self):
        return len(self.sequences)

    def is_soft_core(self, threshold: int) -> bool:
        return len(self.get_genome_names()) >= threshold
    
    def set_soft_core(self, threshold: int):
        self.soft_core = self.is_soft_core(threshold)
        for seq in self.sequences:
            seq.set_soft_core(self.soft_core)
    
    def get_seq_dict(self):
        return {seq.seq_name: seq for seq in self.sequences}
    
    def get_genome_names(self) -> typing.Set[str]:
        return set(seq.get_genome_name() for seq in self.sequences)
    
    def get_sequence_names(self) -> typing.Set[str]:
        return set(seq.seq_name for seq in self.sequences)

    def filter(self, seq_names: list[str]):
        res_coll = copy.deepcopy(self)
        seq_dict = self.get_seq_dict()
        seq_dict = {seq_name: seq_dict[seq_name] for seq_name in seq_names
                if seq_name in seq_dict}
        res_coll.sequences = list(seq_dict.values())
        # return SeqCollection(self.id, list(seq_dict.values()))
        return res_coll

    def align_sequences(self, out_dir = "."):
        '''
        returns an alignment of sequences in synteny block
        '''
        sequences_to_aln = [seq for seq in self.sequences]

        fasta_name = os.path.join(out_dir, "sequences.fasta")
        # SeqIO.write(sequences, fasta_file, "fasta")
        fasta_file = open(fasta_name, "w")
        for seq in sequences_to_aln:
            fasta_file.write(">" + seq.seq_name +"\n")
            fasta_file.write(seq.seq.replace("-", "") + "\n")
            fasta_file.flush()
        fasta_file.close()
        cmdline=["mafft",fasta_name]
        proc = subprocess.run(cmdline,capture_output=True,text=True)
        alignment = proc.stdout
        # print(alignment)
        # os.remove(fasta_name)
        # print(alignment)
        mafft_out = open("mafft.out", "w")
        mafft_out.write(alignment)
        mafft_out.close()
        aligned_seqs = parse_mafft_output("mafft.out")
        # print(aligned_seqs)
        new_seqs = set()
        for seq in sequences_to_aln:
            if seq.seq_name in aligned_seqs:
                seq.seq = aligned_seqs[seq.seq_name]
                new_seqs.add(seq)
            else:
                print("Lack of seq" + seq.seq_name)
        self.sequences = new_seqs

    def to_MAF_block(self, chr_names = None, align = False):
        """
        return a string representing block in MAF format
        """
        if align:
            self.align_sequences()

        seq_lens = np.array([len(seq.seq) for seq in self.sequences])
        # seqs_refound = [seq.refound for seq in self.sequences]
        assert all(seq_lens == seq_lens[0]), f"{self.id},{seq_lens},"
        seq_dict = self.get_seq_dict()
        # if self.aligned:

        maf_str = "a\n"

        if chr_names:
            for maf_seq in sorted(seq_dict.values(), key = lambda x: x.seq_name):
                if maf_seq.seq_name in chr_names:
                    maf_str += maf_seq.MAF_repr() 

        else:
            for maf_seq in sorted(seq_dict.values(), key = lambda x: x.seq_name):
                maf_str += maf_seq.MAF_repr() 

        maf_str += '\n'       

        # else:
        #     sys.exit("MAF block alignment to be implemented")
        
        return maf_str
    
    def convert_coord_system(self, out_format):
        assert out_format in FORMATS_COORD_SYSTEMS
        res_seq_coll = copy.deepcopy(self)
        new_seqs = set()
        for seq in res_seq_coll.sequences:
            # seq.convert_coords(out_format)
            res_seq = seq.convert_coords(out_format)
            # print(res_seq)
            new_seqs.add(res_seq)
        res_seq_coll.sequences = new_seqs
        return res_seq_coll
        # print(self.sequences)
    
    def get_mean_len(self):
        return sum([len(seq) for seq in self.sequences])/len(self)
    
    def get_size(self):
        return len(self)

class Pangenome:
    ### TODO genome lens
    # maybe stores mapping from ids to genome names?
    def __init__(self, seqs_collections: list[SeqCollection]):
        # self.seq_collections = {collection.id: collection for collection in seqs_collections}
        if not all([seq_coll.id for seq_coll in seqs_collections]):
            new_id = 0
            for seq_coll in seqs_collections:
                seq_coll.id = new_id
                new_id += 1
        for seq_coll in seqs_collections:
            clust_id = seq_coll.id
            for seq in seq_coll.sequences:
                seq.cluster_id = clust_id
        self.seq_collections: typing.Set[SeqCollection] = {collection for collection in seqs_collections}
        # self.source_seqs = {seq.name: seq for collection_seqs in seqs_collections for seq in collection_seqs}
        # ??? it would be better to have source seqs assigned to genomes
        # having a set of genomes is necessery for soft_core genome deduction
        # self.genome_names = se t([seq.get_genome_name() for seq_coll in seqs_collections for seq in seq_coll.seq_dict.values()])
        # print(seqs_collections)
        coords_systems = set([seq.coord_system for seq_coll in seqs_collections for seq in seq_coll.sequences])
        assert len(coords_systems) == 1, coords_systems
        self.coord_system = coords_systems.pop()

    def size(self):
        return len(self.seq_collections)

    def get_seq_collections_dict(self):
        return {seq_coll.id: seq_coll for seq_coll in self.seq_collections}
    
    def get_genome_names(self):
        return set([seq.get_genome_name() for seq_coll in self.seq_collections for seq in seq_coll.sequences])
    
    def detect_soft_core(self, threshold: float = 0.95):
        n_genomes_threshold = round(len(self.get_genome_names()) * threshold)
        for seq_coll in self.seq_collections:
            seq_coll.set_soft_core(n_genomes_threshold)

    def get_seq_names(self):
        return set([seq.seq_name for seq_coll in self.seq_collections for seq in seq_coll.sequences])

    def get_all_sequences(self):
        return set([seq for seq_coll in self.seq_collections for seq in seq_coll.sequences])

    def get_sequences_by_seq_name(self):
        res_dict = {seq_name:[] for seq_name in self.get_seq_names()}
        all_seqs = self.get_all_sequences()
        for seq in all_seqs:
            res_dict[seq.seq_name].append(seq)
        
        return res_dict

    def filter(self, seq_names: list[str], preserve_collections=False):
        """
        Results in a model only clusters with chosen sequences are in
        If preserve collection is set to True, whole selected sequence collection
        (together with sequences the model is not filtered by) will be returned.
        """
        seq_names = set(seq_names)
        collection_ids = set()
        filtered_collections = []
        collections_dict = self.get_seq_collections_dict()
        for id, C in collections_dict.items():
            if seq_names.issubset(set(C.seq_dict.keys())):
                collection_ids.add(id)
                if preserve_collections:
                    filtered_collections.append(C)
                else:
                    filtered_collections.append(C.filter(seq_names))
        return Pangenome(filtered_collections)

    def filter_by_genome(self, genome_names: list[str], preserve_collections=False):
        """
        Results in a model only clusters with chosen sequences are in
        If preserve collection is set to True, whole selected sequence collection
        (together with sequences the model is not filtered by) will be returned.
        """
        genome_names = set(genome_names)
        res_model = copy.deepcopy(self)
        collection_ids = set()
        filtered_collections = []
        collections_dict = self.get_seq_collections_dict()
        for id, C in collections_dict.items():
            if genome_names.issubset(set([seq.get_genome_name() for seq in C.sequences])):
                collection_ids.add(id)
                if preserve_collections:
                    filtered_collections.append(C)
                else:
                    seq_names = [seq.seq_name for seq in C.sequences if seq.get_genome_name() in genome_names]
                    filtered_collections.append(C.filter(seq_names))
        res_model.seq_collections = list(filtered_collections)
        # return Pangenome(filtered_collections)
        return res_model

    def get_filtered_vertices_by_strand(self, seq_names, symmetrical_invert=False):
        """"
        vertex - seq collection
        """
        # collections_dict = self.get_seq_collections_dict()
        filtered_collections = {"".join(strand_comb): [] for strand_comb in product(["+","-"], repeat=len(seq_names))}
        for V in self.seq_collections:
            if not set(seq_names).issubset(V.get_sequence_names()):
                continue
            dict_key = ""
            coords = []
            seq_dict = V.get_seq_dict()
            for _seq_name in seq_names:
                seq = seq_dict[_seq_name]
                coords.append((seq.start, seq.end))
                if seq.strand > 0:
                    dict_key += "+"
                else:
                    dict_key += "-"
            filtered_collections[dict_key].append(tuple(coords))
                            
            if symmetrical_invert:
                dict_key = ""
                coords = []
                for seq_name_ in seq_names:
                    seq = seq_dict[seq_name_]
                    coords.append(seq.invert_coords())
                    if seq.strand > 0:
                        dict_key += "-"
                    else:
                        dict_key += "+"
                filtered_collections[dict_key].append(tuple(coords))      

        return filtered_collections

    # def convert_coords_system(self, out_format: str):
    #     assert (out_format in FORMATS_COORD_SYSTEMS)
    #     ### !!!! Mabye dict would be better, as for set attribute you have to complitety update it all,
    #     ### so at once there are two of them
    #     # old_model = copy.deepcopy(self)  
    #     seq_colls = set()
    #     for seq_coll in self.seq_collections:
    #         seq_coll.convert_coord_system(out_format)
    #         res_seq_coll = seq_coll
    #         # print(res_seq_coll)
    #         # print("Res seq")
    #         # for seq in res_seq_coll.sequences:
    #         #     print(seq)
    #         seq_colls.add(res_seq_coll)
    #     self.seq_collections = seq_colls

    #     self.coord_system = out_format

    #     # return old_model   
    def convert_coords_system(self, out_format: str):
        assert (out_format in FORMATS_COORD_SYSTEMS)
        ### !!!! Mabye dict would be better, as for set attribute you have to complitety update it all,
        ### so at once there are two of them
        res_model = copy.deepcopy(self)  
        seq_colls = set()
        for seq_coll in res_model.seq_collections:
            res_seq_coll = seq_coll.convert_coord_system(out_format)
            # res_seq_coll = seq_coll
            # print(res_seq_coll)
            # print("Res seq")
            # for seq in res_seq_coll.sequences:
            #     print(seq)
            seq_colls.add(res_seq_coll)

        res_model.seq_collections = seq_colls
        res_model.coord_system = out_format

        return res_model

    def to_MAF(self, out_file, chr_names=None):
        """
        Writes MAF object into a maf file. Filtering on sequence names
        is possible
        """
        if not self.coord_system == "maf":
            self.convert_coords_system("maf")
        ### alignments have to be prepared if  sequences are not aligned already (check seq lens?)
        out_file = open(out_file,"w")
        out_file.write("##maf version=1 scoring=N/A\n\n")

        if not chr_names:
            for block in self.seq_collections:
                if len(block.sequences) > 0:
                    out_file.write(block.to_MAF_block())
        
        else:
             for block in self.seq_collections:
                if len(block.sequences) > 0:
                    out_file.write(block.to_MAF_block(chr_names))           

        out_file.close()

    def to_GFF(self, out_file):
        if not self.coord_system == "gff":
            self.seq_collections = self.convert_coords_system("gff").seq_collections
            self.coord_system = "gff"
        
        out_file = open(out_file, "w")
        for seq_coll in self.seq_collections:
            coll_id = seq_coll.id
            core_stat = seq_coll.soft_core
            for seq in seq_coll.sequences:
                strand_sign = "+" if seq.strand > 0 else "-"
                out_file.write(f"{seq.seq_name}\tunidentified\tunidentified\t{seq.start}\t{seq.end}\t.\t{strand_sign}\t0\tcluster_id={coll_id};core_status={core_stat}\n")
        
        out_file.close()


    def map_to_gff(self, gff_dir, overlap_threshold = 0.8, sep="."):
        if not self.coord_system == "gff":
            old_coord_cystem = self.coord_system
            self.seq_collections = self.convert_coords_system("gff").seq_collections
            self.coord_system = "gff"
        
        simple_gffs = parse_GFFs_dir(gff_dir)
        genome_cds_coords = simple_gffs.get_genome_cds_coords_dict()
        for seq_coll in self.seq_collections:
            # try:
            #     print(seq_coll.cluster_name)
            # except:
            #     print(seq_coll.id)
            for seq in seq_coll.sequences:
                #sometimes in maf files there will be ancestor sequences that are not
                # present in gff files
                if not seq.seq_name in genome_cds_coords.keys():
                    continue
                for annot in genome_cds_coords[seq.seq_name]:
                    if contains((seq.start, seq.end),(annot.start, annot.end), threshold=overlap_threshold):
                        seq.mapped_annotations.append(annot)
                        # print(seq.annotation_ids)
                        # print(seq.mapped_annotations)
                seq.mapped_annotations = sorted(seq.mapped_annotations, key = lambda x: x.start)
                # print("from_panaroo",seq.annotation_ids)
                # print("mapped",[annot.annotation_id for annot in seq.mapped_annotations])
            # print("-"*40)
        self.seq_collections = self.convert_coords_system(old_coord_cystem).seq_collections
        self.coord_system = old_coord_cystem

    def annotations_to_csv(self, gff_dir, csv_seqs_name = "sequences_summary.csv",csv_annots_name = "annotations_summary.csv", overlap_threshold = 0.7):
        res_csv = open(csv_seqs_name, "w")
        res_csv.write("cluster id,cluster size,core status,seq name,seq start,seq end,seq strand,mapped annotations\n")
        res_csv.flush()
        annots_csv = open(csv_annots_name, "w")
        annots_csv.write("seq name,annotation ID,start,end,strand\n")
        annots_csv.flush()
        self.detect_soft_core()
        self.map_to_gff(gff_dir, overlap_threshold=overlap_threshold)
        for seq_coll in self.seq_collections:
            id = seq_coll.id
            clust_size = len(seq_coll)
            core_status = seq_coll.soft_core
            for seq in seq_coll.sequences:
                # annotation len is also an iimportant aspect, but can be easily retrived from gff file
                # print(seq.seq_name, [ann.annotation_id for ann in seq.mapped_annotations], seq_coll.soft_core)
                annots = ";".join([ann.annotation_id for ann in seq.mapped_annotations])
                res_csv.write(f"{id},{clust_size},{core_status},{seq.seq_name},{seq.start},{seq.end},{seq.strand},{annots}\n")                
                res_csv.flush()
                for annot in seq.mapped_annotations:
                    annots_csv.write(f"{seq.seq_name},{annot.annotation_id},{annot.start},{annot.end},{annot.strand}\n")
                    annots_csv.flush()
        res_csv.close()
        annots_csv.close()

    def get_sequence_adjecency_dict(self):
        """
        For each sequence, which sequence is adjecent to which
        ASSUMES THAT ALL SEQUENCES DO NOT OVERLAP
        """
        adjecency_dict = {}
        # for seq_coll_1 in self.seq_collections:

        #     for seq_coll_2 in self.seq_collections:
        #         if seq_coll_1.id == seq_coll_2.id:
        #             continue
        # model_gff = self.convert_coords_system("gff")
        # seqs_by_seq_name = model_gff.get_sequences_by_seq_name()
        seqs_by_seq_name = self.get_sequences_by_seq_name()
        # print("seqs by seq name", seqs_by_seq_name)
        for seq_name, seqs in seqs_by_seq_name.items():
            # sort sequences by start pos in gff coords - easy to get adjencent fragments
            # from there
            seqs.sort(key = lambda s: s.convert_coords("gff").start)
            # print(seqs)
            # The following works, if sequences are non-overlapping - not the case for original
            # panarpoo input files
            for i in range(0,len(seqs)-1):
                adjecency_dict[seqs[i]] = seqs[i + 1]

            # For the overlapping sewuences case, both overlapiing with start > and end >  should 
            # probably go into adjencent sequences - but for sure not only te fisrt one after the sequence.
            # Nevetheless simple approach will for for 
                # print(seqs[i].convert_coords("gff"), seqs[i+1].convert_coords("gff"))
        # print("adj dict")
        # print(adjecency_dict)     
        return adjecency_dict
    
    def get_cluster_adjency_dict(self):
        sequence_adj_dict = self.get_sequence_adjecency_dict()
        cluster_adj_dict = {seq_coll.id: set() for seq_coll in self.seq_collections}

        for seq, adj_seq in sequence_adj_dict.items():
            # for adj_seq_ in adj_seqs:
            cluster_adj_dict[seq.cluster_id].add(adj_seq.cluster_id)
        return cluster_adj_dict
    
    def get_panaroo_edges(self):
        """
        Returns tuples (seq_coll_1, seq_coll_2) according to Panaroo graph construction:
        'Panaroo builds a full graphical representation of the pangenome, where nodes are
        clusters of orthologous genes (COGs) and two nodes are connected by an edge if they
        are adjacent on a contig in any sample from the population'
        """
        cluster_adj_dict = self.get_cluster_adjency_dict()
        G_edges = []
        for clust_id, adj_clust_ids in cluster_adj_dict.items():
            for adj_id in adj_clust_ids:
                G_edges.append((clust_id, adj_id))
        return G_edges
    
    def get_cluster_graph(self):
        """
        Returns networxx graph constructed according to the snippet from
        the panaroo paper:
        'Panaroo builds a full graphical representation of the pangenome, where nodes are
        clusters of orthologous genes (COGs) and two nodes are connected by an edge if they
        are adjacent on a contig in any sample from the population'
        """
        G_edges = self.get_panaroo_edges()
        return nx.DiGraph(G_edges)

    def get_panaroo_cluster_triplets(self):
        """
        Retrives triplets used by Panaroo to detect structural variants.
        For each genome, for each contig, seq collections are sorted accoring to that genome sequences.
        Blocks have to be paralogs free. For chosen genome, returns list of triplets that 
        are present in this genome
        """
        # filtered_model = self.filter_by_genome([genome])
        # for genome in self.get_genome_names():
        G_clusters = self.get_cluster_graph()
        # cluster_paths = nx.all_simple_paths(G_clusters)
        # print(cluster_paths)
        # cluster_paths = set()
        print("graph constructed")
        cluster_paths = []
        for start_node in G_clusters.nodes:
            for end_node in G_clusters.nodes:
                paths = nx.all_simple_paths(G_clusters, start_node, end_node)
                for path_ in paths:
                    if len(path_) > 2:
                        # cluster_paths.append(path_)
                        for i in range(0, len(path_) - 2):
                            # cluster_paths.add(tuple(path_[i:i+3]))
                            if not tuple(path_[i:i+3]) in cluster_paths:
                                cluster_paths.append(tuple(path_[i:i+3]))

                            # print(cluster_paths)

        return cluster_paths
    
    def triplets_to_csv(self, csv_out):
        triplets = self.get_panaroo_cluster_triplets()
        csv_out = open(csv_out, "w")
        csv_out.write("triplet\n")
        for triplet in triplets:
            csv_out.write(f'{"-".join(map(str, triplet))}\n')
            csv_out.flush()
        csv_out.close()

    def get_panaroo_annotaion_triplets(self, gff_dir):
        self.map_to_gff(gff_dir)
        cluster_triplets = self.get_panaroo_cluster_triplets()
        seq_coll_dict = self.get_seq_collections_dict()

    
    def get_block_sizes_count(self):
        size_counter = {}
        for block in self.seq_collections:
            block_size = block.get_size()
            if block_size in size_counter:
                size_counter[block_size] += 1
            else:
                size_counter[block_size] = 1
        return size_counter
    
    def block_sizes_csv(self, csv_out="blocks_sizes.csv"):
        csv_res = open(csv_out, "w")
        csv_res.write("block size,blocks count\n")
        size_counter = self.get_block_sizes_count()
        for size_, count_ in size_counter.items():
            csv_res.write(f"{size_},{count_}\n")
            csv_res.flush()

        csv_res.close()
    
    def get_block_lens_count(self):
        lens_counter = {}
        for block in self.seq_collections:
            mean_len = round(block.get_mean_len())
            if mean_len in lens_counter:
                lens_counter[mean_len] += 1
            else:
                lens_counter[mean_len] = 1
        return lens_counter

    def block_lens_csv(self, csv_out="blocks_lens.csv"):
        csv_res = open(csv_out, "w")
        csv_res.write("block mean seq len,blocks count\n")
        lens_counter = self.get_block_lens_count()
        for len_, count_ in lens_counter.items():
            csv_res.write(f"{len_},{count_}\n")
            csv_res.flush()

        csv_res.close()
                

    def get_genome_str_presence_absence(self, genome):
        pass
