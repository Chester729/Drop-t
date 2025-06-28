"""
make a global graph (f-graph) from files of fragments and contacts
"""

import sys
import numpy as np
import igraph as ig
import sys
from requirements import make_MboI
import argparse
parser = argparse.ArgumentParser(description="Process some parameters.")


parser.add_argument('-f','--frag_file', required=True, help="frag file")
parser.add_argument('-c','--contact_file', required=True, help="contact file")
parser.add_argument('-o','--graph_file', default = 'combined_bc_contacts_30_frag_weighted.GraphML', help="contact file")
parser.add_argument('-m','--MboI', required=True, help="MboI restriction site file")

args = parser.parse_args()

chrindex={}
for i in range(1,23):
    chrindex[str(i)]=i
chrindex['X']=23
chrindex['Y']=24
def make_MboI(filepath):
    """Load MboI restriction site data"""
    MboI = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            chrom = parts[0]
            if chrom not in chrindex:
                continue
            try:
                MboI[chrom] = np.array([1] + [int(x) for x in parts[1:]])
            except ValueError:
                continue
    return MboI

MboI = make_MboI(args.MboI)

vertex = set()
edge = {}
def make_vertex_edge(bc_frags_file,bc_contacts_file):
    l=0
    vertex = set()
    edge = {}
    with open(bc_frags_file, 'r') as fin:
        for r in fin:
            l += 1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            r = r[1:]
            for i in r:
                ch = i.split('#')[0]
                vertex.add(i)
    with open(bc_contacts_file,'r') as fin:
        for r in fin:
            l+=1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            r = r[1:]
            for i in r:
                frags = i.split(':')
                weight = int(frags[1])
                frags = frags[0]
                frags = frags.split(',')
                frag1 = frags[0][1:]
                frag2 = frags[1][:-1]
                frag1_split = frag1.split('#')
                frag2_split = frag2.split('#')
                ch1 = frag1_split[0]
                ch2 = frag2_split[0]
                loc1 = int(frag1_split[1])
                loc2 = int(frag2_split[1])
                now_pairs = (frag1,frag2)
                if (ch1==ch2 and loc1>loc2) or chrindex[ch1] > chrindex[ch2]:
                    now_pairs = (frag2,frag1)
                try:
                    edge[now_pairs] += 1
                except KeyError:
                    edge[now_pairs] = 1
    print('\n')
    return vertex,edge
vertex,edge = make_vertex_edge(bc_frags_file=args.frag_file,bc_contacts_file=args.contact_file)
print('\n')

G = ig.Graph(directed=False)
G.add_vertices(list(vertex))
G.add_edges(list(edge.keys()),attributes={'weight':list(edge.values())})
del vertex
del edge
G.write_graphml(args.graph_file)
ig.summary(G)
del G