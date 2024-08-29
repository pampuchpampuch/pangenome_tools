from pgtools import panaroo_parser
from pgtools.gff_parser import parse_gff
import os
from Bio.Seq import Seq

panaroo_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset/panaroo_out"
gffs_dir = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset/gff"
panaroo_obj = panaroo_parser.parse_panaroo_output(panaroo_dir, gffs_dir, include_refound=True)

