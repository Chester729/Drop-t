import sys
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
from itertools import combinations
from collections import Counter
from matplotlib.pylab import *
import gzip
import pickle
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['svg.fonttype']='none'
"""
This code calculates the probability of NO common neighbors in the graph for two components randomly sampled from different droplets.
The probability of components from different droplets having common neighbor should be very LOW, so the proportion of HAVING common neighbors reflects the FALSE connection rate for some parameter combinations.
We do this test to determine the parameters of the final algorithm, by testing the probability of no common neighbors for components of different sizes.
"""
import argparse
parser = argparse.ArgumentParser(description="Process some parameters.")
parser.add_argument('-t','--top_frags', required=True,help="read top frags")
parser.add_argument('-m','--MboI', required=True, help="MboI restriction site file")
parser.add_argument('-g','--Graph', required=True, help="Graph file")
parser.add_argument('-c','--component_file',required=True, help="output component file")
parser.add_argument('-d','--droplet_number',default=1000, help="number of droplets for testing",type=int)
parser.add_argument('-sl','--size_low', default=1, help="min component size for testing",type=int)
parser.add_argument('-sh','--size_high', default=100, help="max component size for testing",type=int)
parser.add_argument('-s10','--size_10', default=50, help="test 10 step from this size",type=int)
parser.add_argument('-el','--expand_low', default=0, help="min expand parameter for testing",type=int)
parser.add_argument('-eh','--expand_high', default=10, help="max expand parameter for testing",type=int)
parser.add_argument('-es','--expand_step', default=1, help="step of expand parameter for testing",type=int)
parser.add_argument('-cd','--comp_dict', default=None, help="stored comp_dict file, containing random component pairs")
parser.add_argument('-cdo','--comp_dict_output', default='comp_dict.pkl', help="store sampled comp_dict file")
parser.add_argument('-o','--no_neighbor_prob', default='no_neighbor_prob_dicts.pkl', help="output no neighbor prob dict")
parser.add_argument('--add_chr', action="store_true", help="add chr")
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

def frag1_smaller_frag2(pair):
    frag1_split = pair[0].split('#')
    frag2_split = pair[1].split('#')
    ch1 = frag1_split[0]
    ch2 = frag2_split[0]
    loc1 = int(frag1_split[1])
    loc2 = int(frag2_split[1])
    now_pairs = (pair[0],pair[1])
    if (ch1==ch2 and loc1>loc2) or chrindex[ch1] > chrindex[ch2]:
        now_pairs = (pair[1],pair[0])
    return now_pairs
def expand_component(comp,cut=5):
    comp_set = set(comp)
    for i in comp:
        ch_loc = i.split('#')
        ch = ch_loc[0]
        loc = int(ch_loc[1])
        for j in range(max(0,loc-cut),min(loc+cut+1,len(MboI[ch]))):
            comp_set.add(ch+'#'+str(j))
    return list(comp_set)

#comp_list: length=tuples of lists
with open(args.top_frags,'r') as fin:
    top_frags = fin.readline().split(' ')

def make_ch_frag_dict(comp_list):
    ch_frag_dict = {}
    for i in comp_list:
        i = i.split('#')
        try:
            ch_frag_dict[i[0]].append(int(i[1]))
        except KeyError:
            ch_frag_dict[i[0]] = [int(i[1])]
    for ch in ch_frag_dict.keys():
        ch_frag_dict[ch]  = sorted(ch_frag_dict[ch])
    return ch_frag_dict
"""
Calculate the number of common neighbors for comp1 and comp2
"""
def calculate_common_neighbors(comp_list,G,expand = 0,details = True,top_frags=set()): 
    neighbors_number = []
    neighbor_yes = []
    neighbor_no = []
    neighbor_path = {}
    neighbor_path_weights = {}
    G_name = G.vs['name']
    l=0
    for pairs in comp_list:
        pairs = [expand_component(pairs[0],expand),expand_component(pairs[1],expand)]
        if set(pairs[0]) & set(pairs[1])!=set():
            #print('overlap!',set(pairs[0]),set(pairs[1]))
            l+=1
            continue
        neighbors0 = []
        neighbors1 = []
        for i in pairs[0]:
            try:
                neighbors0+=[G_name[j] for j in G.neighbors(i)]
            except ValueError:
                pass
        for i in pairs[1]:
            try:
                neighbors1+=[G_name[j] for j in G.neighbors(i)]
            except ValueError:
                pass
        neighbors0 = set(neighbors0)-top_frags
        neighbors1 = set(neighbors1)-top_frags
        neighbors_number.append(len(neighbors0&neighbors1))
        if not details:
            l+=1
            continue
        overlaps = neighbors0&neighbors1
        if len(overlaps)>0:
            neighbor_yes.append(pairs)
        else:
            neighbor_no.append(pairs)
        path_set = set()   
        for i in overlaps:      
            frags = set([G_name[j] for j in G.neighbors(i)])
            frags0 = frags&set(pairs[0])  
            frags1 = frags&set(pairs[1])
            for j in frags0:
                path_set.add((j,i,G.es[G.get_eid(j,i)]['weight']))
            for j in frags1:
                path_set.add((i,j,G.es[G.get_eid(i,j)]['weight']))
            for f0 in frags0:
                for f1 in frags1:
                    try:
                        neighbor_path[(l,f0,f1)].add((i,G.es[G.get_eid(f0,i)]['weight'],G.es[G.get_eid(i,f1)]['weight']))
                    except KeyError:
                        neighbor_path[(l,f0,f1)] = {(i,G.es[G.get_eid(f0,i)]['weight'],G.es[G.get_eid(i,f1)]['weight'])}
                    
        weight = 0
        for path in path_set:
            weight+=path[2]
        if weight>0:
            neighbor_path_weights[l] = weight
        l+=1
    if not details:
        return neighbors_number
    return neighbors_number,neighbor_yes,neighbor_no,neighbor_path,neighbor_path_weights

"""
comp_list: list of tuples (comp1,comp2)
"""
def make_comp_list(file,size1,size2,number,start = 3,mode='equal'):
    comp1 = []
    comp2 = []
    l=0
    bc_odd = 0
    now_bc = ''
    with open(file,'r') as fin:
        for r in fin:
            l+=1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            bc = r[0]
            if bc!=now_bc:
                bc_odd = 1-bc_odd
                now_bc = bc
            frags = r[start:]
            frags = list(set(frags)-top_frags)
            if mode=='equal':
                if len(frags)==size1 and bc_odd==1:
                    comp1.append(frags)
                elif len(frags)==size2 and bc_odd==0:
                    comp2.append(frags)
            elif mode=='larger':
                if len(frags)>size1 and bc_odd==1:
                    comp1.append(frags)
                elif len(frags)>size2 and bc_odd==0:
                    comp2.append(frags)
            elif mode=='between':
                if len(frags)>size1[0] and len(frags)<=size1[1] and bc_odd==1:
                    comp1.append(frags)
                elif len(frags)>size2[0] and len(frags)<=size2[1] and bc_odd==0:
                    comp2.append(frags)
            if len(comp1)>=number and len(comp2)>=number:
                break
    comp1 = comp1[:number]
    comp2 = comp2[:number]
    comp_list = [(comp1[i],comp2[i]) for i in range(len(comp1))]
    return comp_list
"""
size_list: [(1,1),(2,2)...(50,50),(51,60),...(91,100)]
Find components of the corresponding sizes
"""
def find_comp_of_size(file,size_list,number,top_frags,special_size):
    l=0
    comp_dict = [{},{}] #odd,even bc
    size_loc_dict = {} #keys: int size, values: locs in size_list
    top_frags = set(top_frags)
    for i in size_list:
        comp_dict[0][i] = []
        comp_dict[1][i] = []
    for i in range(len(size_list)):
        for j in range(size_list[i][0],size_list[i][1]+1):
            size_loc_dict[j] = i
    bc_odd = 0
    now_bc = ''
    if file[-3:]=='.gz':
        fin = gzip.open(file,'rt')
    else:
        fin = open(file,'r')
    for r in fin:
        l+=1
        if l % 10000 == 0:
            sys.stdout.write("\r%d" % l)
            sys.stdout.flush()
        r = r.strip().split(' ')
        bc = r[0]
        if bc!=now_bc:
            bc_odd = 1-bc_odd
            now_bc = bc
        frags = r[1:]
        frags = list(set(frags)-top_frags)
        if len(frags)==0:
            continue
        try:
            now_size = size_list[size_loc_dict[(len(frags))]]
            if len(comp_dict[bc_odd][now_size]) < number:
                comp_dict[bc_odd][now_size].append(frags)
            #print('append!')
        except KeyError:
            if len(frags)<size_list[-1][1]:
                print('Error!\n',len(frags),size_loc_dict)
                return 0
        if l%5000000==0:
            print(list(map(len,comp_dict[0].values())))
            print(list(map(len,comp_dict[1].values())))
        if sum([len(comp_dict[0][size])>=number for size in special_size] + [len(comp_dict[1][size])>=number for size in special_size])>=2*len(special_size):
            break
    return comp_dict

def make_comp_list_from_dict(comp_dict,size1,size2,number):
    try:
        return [(comp_dict[0][size1][i],comp_dict[1][size2][i]) for i in range(number)]
    except IndexError:
        number = min(len(comp_dict[0][size1]),len(comp_dict[1][size2]))
        #print(number)
        return [(comp_dict[0][size1][i],comp_dict[1][size2][i]) for i in range(number)]
def calculate_no_neighbor_prob(comp_dict,size1,size2,expand_list,G,number=1000,details=True,top_frags=set()): 
    no_neighbor_prob_big = {}
    for i in range(len(size1)):
        sys.stdout.write("\r%d" % i)
        sys.stdout.flush()
        comp_list = make_comp_list_from_dict(comp_dict,size1[i],size2[i],number)
        for j in expand_list:
            common_neighbors = calculate_common_neighbors(comp_list,G,expand = j,details=details,top_frags=top_frags)
            try:
                no_neighbor_prob_big[j].append(sum(np.array(common_neighbors)==0))
            except KeyError:
                no_neighbor_prob_big[j] = [sum(np.array(common_neighbors)==0)]
    return no_neighbor_prob_big


import igraph as ig
Gf = ig.Graph.Read_GraphML(args.Graph)

#comp_file = sys.argv[2]
if args.comp_dict==None:
    comp_dict = find_comp_of_size(args.component_file,[(i,i) for i in range(1,args.size_10+1)]+[(i*10+1,i*10+10) for i in range(args.size_10//10,args.size_high//10)],args.droplet_number,top_frags,[(args.size_10,args.size_10),(args.size_high-9,args.size_high),(args.size_10-1,args.size_10-1)])
    with open(args.comp_dict_output,'wb') as f:
        pickle.dump(comp_dict,f)
else:
    with open(args.comp_dict,'rb') as f:
        comp_dict = pickle.load(f)
expand_list = [i for i in range(args.expand_low,args.expand_high+1,args.expand_step)]

no_neighbor_prob = calculate_no_neighbor_prob(comp_dict,[(i,i) for i in range(1,args.size_10+1)],[(i,i) for i in range(1,args.size_10+1)],expand_list,Gf,details=False,top_frags=set(top_frags),number = args.droplet_number)
no_neighbor_prob_1 = calculate_no_neighbor_prob(comp_dict,[(1,1) for i in range(1,args.size_10+1)],[(i,i) for i in range(1,args.size_10+1)],expand_list,Gf,details=False,top_frags=set(top_frags),number = args.droplet_number)
no_neighbor_prob_5 = calculate_no_neighbor_prob(comp_dict,[(5,5) for i in range(1,args.size_10+1)],[(i,i) for i in range(1,args.size_10+1)],expand_list,Gf,details=False,top_frags=set(top_frags),number = args.droplet_number)


no_neighbor_prob_big = calculate_no_neighbor_prob(comp_dict,[(i*10+1,i*10+10) for i in range(args.size_10//10,args.size_high//10)],[(i*10+1,i*10+10) for i in range(args.size_10//10,args.size_high//10)],expand_list,Gf,details=False,top_frags=set(top_frags),number = args.droplet_number)
no_neighbor_prob_big_1 = calculate_no_neighbor_prob(comp_dict,[(1,1) for i in range(args.size_10//10,args.size_high//10)],[(i*10+1,i*10+10) for i in range(args.size_10//10,args.size_high//10)],expand_list,Gf,details=False,top_frags=set(top_frags),number = args.droplet_number)
no_neighbor_prob_big_5 = calculate_no_neighbor_prob(comp_dict,[(5,5) for i in range(args.size_10//10,args.size_high//10)],[(i*10+1,i*10+10) for i in range(args.size_10//10,args.size_high//10)],expand_list,Gf,details=False,top_frags=set(top_frags),number = args.droplet_number)
import pickle
with open(args.no_neighbor_prob,'wb') as f:
    pickle.dump([no_neighbor_prob,no_neighbor_prob_1,no_neighbor_prob_5,no_neighbor_prob_big,no_neighbor_prob_big_1,no_neighbor_prob_big_5],f)