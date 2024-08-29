import typing
from itertools import product
import numpy as np

FORMATS_COORD_SYSTEMS = {
    "maf": {"coord_system": ("0-based", "half-open"), "rev_strand_coords": "rev_strand"},
    "gff": {"coord_system": ("1-based", "half-open"), "rev_strand_coords": "forward_strand"},
    "gfa": {"coord_system": ("0-based", "half-open"), "rev_strand_coords": "forward_strand"},
    # is panaroo coord system the same for refound and regular (gff) genes - NO, refound is 0-based?
    "panaroo": {"coord_system": ("1-based", "half-open"), "rev_strand_coords": "forward_strand"},
    "panaroo_refound":  {"coord_system": ("0-based", "half-open"), "rev_strand_coords": "forward_strand"}
    }

def convert_coords(start: int, end: int, strand: int, src_size: int, in_format: str, out_format: str = "maf"):
    ### !!! currently only conversion to MAF format coord system is implemented.
    assert out_format == "maf"
    coord_system = FORMATS_COORD_SYSTEMS[in_format]["coord_system"]

    res_start = start if coord_system[0] == "0-based" else start - 1
    res_end = end if coord_system[1] == "half-open" else end + 1

    if strand > 0 or FORMATS_COORD_SYSTEMS[in_format] == FORMATS_COORD_SYSTEMS[out_format]:
        return res_start, res_end        
    
    rev_res_start = src_size - res_end
    rev_res_end = src_size - res_start

    return rev_res_start, rev_res_end
        
# based i closed sprowadza sie tylko do tego, że przysplicowaniu w pythonie:
# dla + strand [start - based, end +1 if fully-closed]
# wtedy długość to ...
# dla - strand 

# class SourceSeq:
#     """
#     Represents a bigger Source fragment the sequences (BaseSeq) are positioned on.
#     Usually a chromosome or a contig
#     """
#     def __init__(self, id: int, seq_name: str, size: int, seq: str = ""):
#         self.id = id
#         self.name = seq_name
#         self.size = size
#         self.seq = seq
    
#     def get_genome_name(self, sep="."):
#         return self.name.split(".")[0]
    
class BaseSeq:
    """
    Represent a genomic sequence as a fragment of bigger Source fragment
    (usually chromosome or a contig)
    As the MAF files were commonly used through out the analysis, default coordination system
    for this class is the one used in MAF files (zero-based, half open and coords for "-" strand
    are given in respect to the reversed strand)
    """
    def __init__(self, seq_name: str, start: int , end: int, strand: int, src_size: int, in_format, seq: str = None, coord_system="maf"):
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
        self.in_core: bool = None
        
        ### Add coord system?
        ### For which strands are coords given in the case of - strand

    def invert_coords(self):
        return self.src_size - self.end, self.src_size - self.start

    def get_genome_name(self, sep="."):
        return self.seq_name.split(sep)[0]
    
    def get_contig_name(self, sep = "."):
        return self.seq_name.split(sep)[1]

    def set_core(self, is_core: bool):
        self.in_core = is_core      
    
    def __print__(self):
        return (f"{self.seq_name}: {self.strand} - {self.end}, {self.strand}")

    def __len__(self):
        """
        Assumes 0-based half-open coord system
        """
        return self.end - self.start

class SeqCollection:
    def __init__(self, id: int, seq_list: list[BaseSeq]):
        self.id: int = id
        # seq_dict = {seq.seq_name: seq for seq in seq_list}
        seq_dict = {seq for seq in seq_list}        
        self.sequences: typing.Set[BaseSeq] = set(seq_list)
        # self.genome_names: typing.Set[str] = set(seq.get_genome_name() for seq in seq_dict.values())
        self.core: bool = None

    def is_core(self, threshold: int) -> bool:
        return len(self.genome_names) >= threshold
    
    def set_core(self, threshold: int):
        self.core = self.is_core(threshold)
        for seq in self.sequences:
            seq.set_core(self.core)
    
    def get_seq_dict(self):
        return {seq.seq_name: seq for seq in self.sequences}
    
    def get_genome_names(self) -> typing.Set[str]:
        return set(seq.get_genome_name() for seq in self.sequences)
    
    def get_sequence_names(self) -> typing.Set[str]:
        return set(seq.seq_name for seq in self.sequences)

    def filter(self, seq_names: list[str]):
        seq_dict = self.get_seq_dict()
        seq_dict = {seq_name: seq_dict[seq_name] for seq_name in seq_names
                if seq_name in seq_dict}
        return SeqCollection(self.id, list(seq_dict.values()))

    def to_MAF(self):
        seq_lens = np.array([len(seq) for seq in self.sequences])
        assert all(seq_lens == seq_lens[0])
    
    def get_alignment(self):
        pass

class Pangenome:
    ### TODO genome lens
    # maybe stores mapping from ids to genome names?
    def __init__(self, seqs_collections: list[SeqCollection]):
        # self.seq_collections = {collection.id: collection for collection in seqs_collections}
        self.seq_collections: typing.Set[SeqCollection] = {collection for collection in seqs_collections}
        # self.source_seqs = {seq.name: seq for collection_seqs in seqs_collections for seq in collection_seqs}
        # ??? it would be better to have source seqs assigned to genomes
        # having a set of genomes is necessery for core genome deduction
        # self.genome_names = se t([seq.get_genome_name() for seq_coll in seqs_collections for seq in seq_coll.seq_dict.values()])
        coords_systems = set([seq.coord_system for seq_coll in seqs_collections for seq in seq_coll.sequences])
        assert len(coords_systems) == 1
        self.coord_system = coords_systems.pop()

    def get_seq_collections_dict(self):
        return {seq_coll.id: seq_coll for seq_coll in self.seq_collections}
    
    def get_genome_names(self):
        return set([seq.get_genome_name() for seq_coll in self.seq_collections for seq in seq_coll.sequences])
    
    def detect_core(self, threshold: float = 0.95):
        n_genomes_threshold = round(len(self.genome_names) * threshold)
        for seq_coll in self.seq_collections:
            seq_coll.set_core(n_genomes_threshold)

    def filter(self, seq_names: list[str]):
        seq_names = set(seq_names)
        collection_ids = set()
        filtered_collections = []
        collections_dict = self.get_seq_collections_dict()
        for id, C in collection_ids.items():
            if seq_names.issubset(set(C.seq_dict.keys())):
                collection_ids.add(id)
                filtered_collections.append(C.filter(seq_names))
        return Pangenome(filtered_collections)

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
    
    def convert_coords_system(self, in_format: str, out_format: str):
        assert in_format and out_format in FORMATS_COORD_SYSTEMS
        pass

    def to_MAF(self):
        if self.coord_system == "maf":
            pass
        else:
            #conversion to maf system coords
            pass