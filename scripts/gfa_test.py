from pgtools.gfa_parser import parse_gfa1
import pickle
gfa = "/home/pampuch/studia/magisterka/test_data/klebsiella_subset_old/klebsiella.gfa"
gfa = parse_gfa1(gfa)

gfa_dump = open("gfa.pickle", "w")
pickle.dump(gfa, gfa_dump)
gfa_dump.close()

gfa_loaded = pickle.load("gfa.pickle")

for v in gfa_loaded.seq_collections():
    print(v.id)