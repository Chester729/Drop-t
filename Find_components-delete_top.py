import sys
import numpy as np
import igraph as ig
import argparse
"""
Find connected components, directly using induced subgraph of each droplet
We do this for further testing parameters by "calculate_no_neighbor_prob.py".
"""
parser = argparse.ArgumentParser(description="Process some parameters.")
parser.add_argument('-g','--Graph', required=True, help="Graph file")
parser.add_argument('-t','--top_frag_dict',required=True, help="read top frag dict")
parser.add_argument('-f','--frag_file', required=True,help="frag file")
parser.add_argument('-c','--component_file',default = 'Components_total-delete_top', help="output component file")
args = parser.parse_args()


G = ig.Graph.Read_GraphML(args.Graph)



top_frag_dict = np.load(sys.argv[4],allow_pickle=True).item()

fout = open(args.component_file,'w')
with open(args.frag_file,'r') as fin:
    l=0
    for r in fin:
        l+=1
        if l % 1000 == 0:
            sys.stdout.write("\r%d" % l)
            sys.stdout.flush()
        bc_frags = {}
        r = r.strip().split(' ')
        bc = r[0]
        r = r[1:]
        bc_frags[bc] = []
        for i in r:
            bc_frags[bc].append(i)
        G_sub = G.induced_subgraph(bc_frags[bc])
        G_sub.delete_vertices(list(top_frag_dict[bc].keys()))
        G_sub.add_vertices(list(top_frag_dict[bc].keys()))         #delete the "top frags"
        for i in top_frag_dict[bc].keys():
            G_sub.add_edges(top_frag_dict[bc][i])
        subname = G_sub.vs['name']
        components_test = G_sub.components()
        for comp in components_test:
            fout.write(bc)
            G_sub_sub = G_sub.induced_subgraph(comp)
            for i in comp:
                fout.write(' ' + subname[i])
            fout.write('\n')
    fout.close()