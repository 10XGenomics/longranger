import ctypes
from ctypes import cdll
import os

lib = cdll.LoadLibrary("target/release/libtra.so")

class RustGraph(ctypes.Structure):
    pass

graph_p = ctypes.POINTER(RustGraph)

lib.make_graph.restype = graph_p
lib.write_graph.argtypes = [graph_p]
lib.read_graph.restype = graph_p
lib.num_edges.argtypes = [graph_p]
lib.delete_graph.argtypes = [graph_p]
lib.create_indexed_graph.restype = graph_p
lib.create_indexed_graph.argtypes = [graph_p, ctypes.c_char_p]
lib.graph_contains_seq.restype = ctypes.c_bool
lib.graph_contains_seq.argtypes = [graph_p, ctypes.c_char_p]

lib.find_graph_edges.restype = ctypes.c_char_p
lib.find_graph_edges.argtypes = [graph_p, ctypes.c_char_p]

g = lib.make_graph("./test_data/test_small_graph")
lib.write_graph(g, "test.bin.gz")
new_g = lib.read_graph("test.bin.gz")
print 'Edges before and after', lib.num_edges(g), lib.num_edges(new_g)

files = ['b' + str(i) + '.bin.gz' for i in range(64)]

for i in range(64):
    lib.compute_bwt_bucket(new_g, i, files[i])

lib.merge_bwt_buckets((ctypes.c_char_p * len(files))(*files), len(files), 2, "bin.gz")
indexed_g = lib.create_indexed_graph(new_g, 'bin.gz')
print lib.num_edges(indexed_g)

for file in files:
    os.remove(file)

lib.write_graph(indexed_g, 'test.bin.gz')
new_g = lib.read_graph('test.bin.gz')

print 'test some sequences'
seqs = ['AAATTAAAGAAGTTTTCC', "CAGTTCCATATGGCT", 'CCCCCCCCCCC', 'AAA']

for seq in seqs: 
    print seq, lib.graph_contains_seq(new_g, seq)

for seq in seqs:
    res = str(lib.find_graph_edges(new_g, seq))
    print seq, res.split(',')

lib.delete_graph(indexed_g)
lib.delete_graph(new_g)
lib.delete_graph(g)

