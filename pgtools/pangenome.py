import typing
from itertools import product
config_formats = {
    "maf": {"coord_system": ("0-based", "half-open"), "rev_strand_coords": "rev_strand"},
    "gff": {"coord_system": ("1-based", "fully-closed"), "rev_strand_coords": "forward_strand"},
    "gfa": {None: None},
    "panaroo": {None: None}
}

# based i closed sprowadza sie tylko do tego, że przyspicowaniu w pythonie:
# dla + strand [start - based, end +1 if fully-closed]
# wtedy długość to ...
# dla - strand 


class BaseSeq:
    """
    Represent a genomic sequence as a fragment of bigger structural fragment
    (usually chromosome or a contig)
    """
    def __init__(self, seq_name: str, start: int , end: int, strand: int):
        self.name = seq_name
        self.start = start
        self.end = end
        assert strand == -1 or strand == 1
        self.strand = strand

        ### Add coord system?
        ### For which strands are coords given in the case of - strand
    
    def convert_coords():
        pass
    
    def __print__(self):
        return (f"{self.name}: {self.strand} - {self.end}, {self.strand}")

    def __len__(self):
        """
        Assumes 0-based half-open coord system
        """
        return self.end - self.start

class SeqCollection:
    def __init__(self, id: int, seq_dict:  typing.Dict[str, str]):
        self.id = id
        self.seq_dict = seq_dict
    
    def filter(self, seq_names: list[str]):
        seq_dict = {seq_name: self.seq_dict[seq_name] for seq_name in seq_names
                if seq_name in self.seq_dict.keys()}
        return SeqCollection(self.id, seq_dict)       

class Pangenome:
    ### TODO genome lens
    # maybe stores mapping from ids to genome names?
    def __init__(self, seqs_collections: list[SeqCollection]):
        self.seq_collections = {collection.id: collection for collection in seqs_collections}

    def filter(self, seq_names: list[str]):
        seq_names = set(seq_names)
        collection_ids = set()
        filtered_collections = []
        for id, C in self.seq_collections.items():
            if seq_names.issubset(set(C.seq_dict.keys())):
                collection_ids.add(id)
                filtered_collections.append(C.filter(seq_names))
        return Pangenome(filtered_collections)

    def get_filtered_vertices_by_strand(self,seq_names):
        filtered_collections = {"".join(strand_comb): [] for strand_comb in product(["+","-"], repeat=len(seq_names))}
        for V in self.seq_collections.values():
            if not set(seq_names).issubset(set(V.seq_dict.keys())):
                continue
            dict_key = ""
            coords = []
            for _seq_name in seq_names:
                seq = V.seq_dict[_seq_name]
                coords.append((seq.start, seq.end))
                if seq.strand > 0:
                    dict_key += "+"
                else:
                    dict_key += "-"

            filtered_collections[dict_key].append(tuple(coords))
        
        return filtered_collections