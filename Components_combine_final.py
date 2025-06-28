import sys
import numpy as np
import igraph as ig
import math
import bisect
import sys
import numpy as np
import igraph as ig
from itertools import combinations
from collections import Counter


import argparse
parser = argparse.ArgumentParser(description="Process some parameters.")
parser.add_argument('-m','--MboI', required=True, help="MboI restriction site file")
parser.add_argument('-g','--Graph', required=True, help="Graph file")
parser.add_argument('-t','--top_frag_dict', required=True,help="top_frag_dict.npy")
parser.add_argument('-d','--droplets',required=True,help='fragments in each droplet')
parser.add_argument('-x', type=int, default=50, help="max size threshold for connected components")
parser.add_argument('-y', type=int, default=8, help="max single component size threshold") 
parser.add_argument('-z', type=int, default=5, help="expansion cutoff for neighborhood")
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
G = ig.Graph.Read_GraphML(args.Graph)
top_frag_dict = np.load(args.top_frag_dict,allow_pickle=True).item()
def get_neighbor_weights(graph,node,G_name,cut=5):
    try:
        neighbors = graph.neighbors(node,mode=ig.ALL)
    except ValueError:
        return set(),{}
    neighbors = [G_name[i] for i in neighbors]
    neighbor_weights = [graph.es[graph.get_eid(node,neighbor)]['weight'] for neighbor in neighbors]
    weight_dict = dict([(neighbors[i],neighbor_weights[i]) for i in range(len(neighbors))])
    return set(neighbors),weight_dict
def expand_component(comp,cut=args.z):
    comp_set = set(comp)
    for i in comp:
        ch_loc = i.split('#')
        ch = ch_loc[0]
        loc = int(ch_loc[1])
        for j in range(max(0,loc-cut),min(loc+cut+1,len(MboI[ch]))):
            comp_set.add(ch+'#'+str(j))
    return list(comp_set)
def find_components_and_expand(r,G,expand):
    bc_frags_bc = []
    r = r.strip().split(' ')
    bc = r[0]
    r = r[1:]
    for i in r:
        bc_frags_bc.append(i)
    bc_frags_set = set(bc_frags_bc)
    G_sub = G.induced_subgraph(bc_frags_bc)
    G_sub.delete_vertices(list(top_frag_dict[bc].keys()))
    subname = G_sub.vs['name']
    components0 = G_sub.components()
    components_new = list(components0)
    components0 = [[subname[i] for i in comp] for comp in components0]
    components0_expand = [expand_component(i,expand) for i in components0]        #扩充了components，每个都变成周围5个
    return components0,components0_expand,components_new,bc,subname
def get_neighbors_and_weights(components0,components0_expand,G,G_name,expand):
                #找components的所有neighbor
        
    comp_node_neighbor_dict = {} #keys: components ,values: dict{keys:node,values:[set(neighbors),weight_dict]}
    comp_size_dict = {}  #keys: size, values:all components of the same size
    comp_neighbor_dict = {} #keys: components, values: neighbor set
    all_neighbor_weights_dict = {}   #keys:neighbors, value:dict{keys:comp,values:all weights of this neighbor from this comp}
    
    for i in range(len(components0)):
        try:
            comp_size_dict[len(components0[i])].append(i)
        except KeyError:
            comp_size_dict[len(components0[i])] = [i]
        comp_node_neighbor_dict[i] = dict([(node,get_neighbor_weights(G,node,G_name,cut = expand)) for node in components0_expand[i]])
        comp_neighbor_dict[i] = set()
        for j in comp_node_neighbor_dict[i].keys():
            comp_neighbor_dict[i] = comp_neighbor_dict[i]|comp_node_neighbor_dict[i][j][0]   #combine all neighbor set

    for comp in comp_node_neighbor_dict.keys():
        for node in comp_node_neighbor_dict[comp].keys():
            for neighbor in comp_node_neighbor_dict[comp][node][1].keys():
                try:
                    all_neighbor_weights_dict[neighbor][comp]+=comp_node_neighbor_dict[comp][node][1][neighbor]
                except KeyError:
                    try:
                        all_neighbor_weights_dict[neighbor][comp]=comp_node_neighbor_dict[comp][node][1][neighbor]
                    except KeyError:
                        all_neighbor_weights_dict[neighbor] = {}
                        all_neighbor_weights_dict[neighbor][comp]=comp_node_neighbor_dict[comp][node][1][neighbor]
    return comp_size_dict,comp_neighbor_dict,all_neighbor_weights_dict
def connection(components_new,comp_size_dict,comp_neighbor_dict,all_neighbor_weights_dict,max_size,max_single):
    sorted_size = sorted(list(comp_size_dict.keys()))
    ii=0
    while ii<len(sorted_size):
        i = sorted_size[ii]    #from size=1
        for comp_num in range(len(comp_size_dict[i])):  #get all comps of size=i
            comp = comp_size_dict[i][comp_num]
            if comp==-1:
                continue
            comp_overlap_dict = {} #record all overlap neighbors
            larger_components = []
            for size in sorted_size:
                if size<i:
                    continue
                for comp2 in comp_size_dict[size]:  #for all components
                    if comp2==comp or comp2==-1:
                        continue
                    overlap_neighbors = comp_neighbor_dict[comp]&comp_neighbor_dict[comp2]
                    if overlap_neighbors!=set():
                        comp_overlap_dict[comp2] = overlap_neighbors
            if comp_overlap_dict=={}:
                continue
            overlap_len = [len(comp_overlap_dict[x]) for x in comp_overlap_dict.keys()]
            candidate_comps = [x for x in comp_overlap_dict.keys() if len(comp_overlap_dict[x]) == max(overlap_len)] #with most neighbor
            candidate_size = np.array([len(components_new[c]) for c in candidate_comps])
            candidate_comps = list(np.array(candidate_comps)[np.where(candidate_size == np.min(candidate_size))[0]]) #the smallest one
            candidate_comps = [now_comp for now_comp in candidate_comps if components_new[now_comp]!=[]]
            sum_weights = np.array([sum([all_neighbor_weights_dict[neighbor][now_comp] + all_neighbor_weights_dict[neighbor][comp] for neighbor in comp_overlap_dict[now_comp]]) for now_comp in candidate_comps])
            sum_weights = sum_weights/np.sum(sum_weights)
            picked_comp = np.random.choice(candidate_comps,p=sum_weights)
            new_length = len(components_new[picked_comp])+len(components_new[comp])
            if new_length>max_size or (len(components_new[picked_comp])>max_single and len(components_new[comp])>max_single):
                continue
#start connection
            comp_size_dict[len(components_new[picked_comp])][comp_size_dict[len(components_new[picked_comp])].index(picked_comp)]=-1  #change the dict for lens
            comp_size_dict[i][comp_num] = -1

            try:
                comp_size_dict[new_length].append(picked_comp)
            except KeyError:
                comp_size_dict[new_length] = [picked_comp]
                bisect.insort_left(sorted_size, new_length)
            components_new[picked_comp] += components_new[comp] #connection, change the component
            comp_neighbor_dict[picked_comp] = comp_neighbor_dict[picked_comp]|comp_neighbor_dict[comp] #change the dict of neighbors
            components_new[comp] = []
            for neighbor in comp_neighbor_dict[comp]:
                try:
                    all_neighbor_weights_dict[neighbor][picked_comp] += all_neighbor_weights_dict[neighbor][comp]
                except KeyError:
                    all_neighbor_weights_dict[neighbor][picked_comp] = all_neighbor_weights_dict[neighbor][comp]
                del all_neighbor_weights_dict[neighbor][comp]
            del comp_neighbor_dict[comp]
        ii+=1
    return components_new
def write_dLHCC(bc,components_new,fout,subname):
    for comp in components_new:
        if comp==[]:
            continue
        comp_add = [subname[i] for i in comp]
        fout.write(bc)
        for i in comp_add:
            fout.write(' ' + i)
        fout.write('\n')
def process_single_droplet(r,G,top_frag_dict,G_name,fout,max_size,max_single,expand):
    components0,components0_expand,components_new,bc,subname = find_components_and_expand(r,G,expand)
    comp_size_dict,comp_neighbor_dict,all_neighbor_weights_dict = get_neighbors_and_weights(components0,components0_expand,G,G_name,expand)
    components_new = connection(components_new,comp_size_dict,comp_neighbor_dict,all_neighbor_weights_dict,max_size,max_single)
    write_dLHCC(bc,components_new,fout,subname)
def process_all_droplets(args):
    file_suffix = args.droplets.split('_')[-1]
    output_dir = '/'.join(args.droplets.split('/')[:-1])
    fout_name = output_dir+'/Components_combine_local_expand_'+str(args.x)+'_'+str(args.y)+'_'+str(args.z)+'_'+file_suffix
    G_name = G.vs['name']
    max_size = int(args.x)      #The max size allowed for the ligation products
    max_single = int(args.y)    #The max size allowed for each component
    expand = int(args.z)                                    #expand 
    with open(args.droplets,'r') as fin,open(fout_name,'w') as fout:
        l=0
        for r in fin:
            l += 1
        #prepare，make induced subgraph, find components
            if l % 50 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            process_single_droplet(r,G,top_frag_dict,G_name,fout,max_size,max_single,expand)
process_all_droplets(args)