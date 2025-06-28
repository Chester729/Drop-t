import numpy as np
import pybedtools
import re
class fragment_analysis:
    def fragment_length(frag,MboI):
        frag = frag.split('#')
        return MboI[frag[0]][int(frag[1])+1] - MboI[frag[0]][int(frag[1])]
    def get_fragment_sequence(frag,genome='hg19'):
        frag = frag.split('#')
        if genome=='hg19':
            genome = '/home/Reference_Genome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
        a = pybedtools.BedTool(f"chr{frag[0]} {MboI[frag[0]][int(frag[1])]} {MboI[frag[0]][int(frag[1])+1]}", from_string=True)
        fasta = pybedtools.example_filename(genome)
        a = a.sequence(fi=fasta)
        print(open(a.seqfn).read())
class complex_analysis:
    def sort_frags(frags):
        frags = list(frags)
        frags = [i.split('#') for i in frags]
        frags = sorted(frags,key=lambda x:(chrindex[x[0]],int(x[1])))
class merged_nodups_analysis:
    def get_matches_cigar(cigar,chrs='[SMHID]',reverse=False):
        nums = [int(x) for x in re.split(chrs,cigar) if x]
        chars = re.findall(chrs,cigar)
        cigar_dict = {}
        m = 0
        s=0
        matching=False
        if reverse:
            chars.reverse()
            nums.reverse()
        mismatches = 0
        mismatch_words = 'IDX'
        for i in range(len(chars)):
            if chars[i] in mismatch_words:
                mismatches+=nums[i]
            if chars[i]=='M':
                m += nums[i]
                matching=True
            elif not matching:
                s += nums[i]
        return m,s,mismatches
class math:
    def cv(x):
        return np.std(x)/np.mean(x)
    def filt_outliers(x,nstd = 5):
        x_mean = np.mean(x)
        x_std = np.std(x)
        x_min = x_mean-nstd*x_std
        x_max = x_mean+nstd*x_std
        x[x<x_min] = x_min
        x[x>x_max] = x_max
        return x
    def mean_std_normalize(x):
        x_mean = np.mean(x)
        x_std  = np.std(x)
        if np.std(x)==0:
            return x
        return (x-x_mean)/x_std
    def max_min_normalize(x):
        xmax = np.max(x)
        xmin = np.min(x)
        if xmax==xmin:
            return x
        return [(i-xmin)/(xmax-xmin) for i in x]
class dict_analysis:
    def counter_dict2ratio(counter_dict):
        new_dict = {}
        counter_dict = dict(counter_dict)
        sums = sum(counter_dict.values())
        for i in counter_dict.keys():
            new_dict[i] = counter_dict[i]/sums
        return new_dict
class set_analysis:
    def Jaccard(A,B):
        C = A & B
        return len(C)/(len(A)+len(B)-len(C))