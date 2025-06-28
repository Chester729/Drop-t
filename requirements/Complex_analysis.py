import gzip
import numpy as np
import sys
from .make_MboI import make_MboI
class complex_analysis:
    def __init__(self, file_name, genome='hg19',get_MboI=False):
        self.file_name = file_name
        self.genome = genome
        if self.file_name.endswith('.gz'):
            self.fin = gzip.open(self.file_name, 'rt')
        else:
            self.fin = open(self.file_name, 'r')
        
        # Read the first line and determine the start index
        line1 = self.fin.readline().strip().split(' ')
        self.start = 0
        for i in range(len(line1)):
            if '#' not in line1[i]:
                self.start += 1
            else:
                break
        
        # Determine if we need to add 'chr'
        self.add_chr = line1[-1][:2].lower() == 'ch'
        
        # Initialize MboI
        #print(self.genome,self.add_chr,line1[-1][:2].lower())
        if get_MboI:
            self.chrindex,self.MboI = make_MboI(self.genome, self.add_chr)
        self.fin.close()
    def set_MboI(self,MboI):
        self.MboI = MboI
    def set_chrindex(self,chrindex):
        self.chrindex = self.chrindex
    def reset_file(self):
        if self.file_name.endswith('.gz'):
            self.fin = gzip.open(self.file_name, 'rt')
        else:
            self.fin = open(self.file_name, 'r')
    def get_complex_num(self):
        self.reset_file()
        l = 0
        for r in self.fin:
            l += 1
        self.complex_num = l
    def set_start(self,start):
        self.start = start
    def get_fragnums(self):
        self.reset_file()
        l = 0
        self.fragnums = []
        for r in self.fin:
            l += 1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')[self.start:]
            self.fragnums.append(len(r[1:]))
        self.fragnums = np.array(self.fragnums)
    
    def get_sum_frag(self, fragnum_cut=1, size_cut=1000):
        self.reset_file()
        self.sum_frag = []
        l = 0
        for r in self.fin:
            l += 1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            s = 0
            r = r.strip().split(' ')[self.start:]
            if len(r) <= fragnum_cut:
                continue
            ch = [i.split('#')[0] for i in r]
            try:
                loc = [i.split('#')[1] for i in r]
            except IndexError:
                print(l, r)
                continue
            if len(ch) != len(loc):
                print(self.file_name, 'Error!\n')
            for i in range(len(ch)):
                try:
                    xx = self.MboI[ch[i]][int(loc[i]) + 1] - self.MboI[ch[i]][int(loc[i])]
                except (KeyError, IndexError, ValueError):
                    continue
                s += xx
                if xx < 0:
                    print(ch[i], loc[i], r)
                    return 0
            if s > size_cut:
                self.sum_frag.append(s)
        self.fin.close()
        self.sum_frag = np.array(self.sum_frag)
    def split_by_chr(self,query_chr,ffout,mode = 'all',ffout_other=None,ffout_hybrid=None):
        query_chr = set(query_chr)
        self.reset_file()
        fout = open(ffout,'w')
        if ffout_other!=None:
            fout_other = open(ffout_other,'w')
        if ffout_hybrid!=None:
            fout_hybrid = open(ffout_hybrid,'w')
        l=0
        for r in self.fin:
            l += 1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r_org = r
            r = r.strip().split(' ')[self.start:]
            ch = set([i.split('#')[0] for i in r])
            if mode=='overlap':  #只要含有这个物种就行
                if ch&query_chr!=set():
                    fout.write(r_org)
            elif mode=='all':   #全是这个物种
                if ch-query_chr ==set():
                    fout.write(r_org)
            elif mode=='separate':
                if ch-query_chr ==set():
                    fout.write(r_org)
                elif ch&query_chr==set():
                    fout_other.write(r_org)
                else:
                    fout_hybrid.write(r_org)
    def make_num_bc_frag(self):
        num_bc_frag = {}
        l=0
        self.reset_file()
        for r in self.fin:
            l += 1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            bc = r[0]
            for f in r[self.start:]:
                if f not in num_bc_frag.keys():
                    num_bc_frag[f]=0
        print('\n',len(num_bc_frag))
        l=0
        self.reset_file()
        for r in self.fin:
            l += 1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            bc = r[0]
            for f in r[self.start:]:
                try:
                    num_bc_frag[f]+=1
                except KeyError:
                    pass
        self.num_bc_frag = num_bc_frag
    def find_lines_of_region(self,ffout,region,is_contacts = False):
        l=0
        with open(ffout,'w') as fout:
            for r in self.fin:
                l += 1
                if l % 10000 == 0:
                    sys.stdout.write("\r%d" % l)
                    sys.stdout.flush()
                r_org = r
                r = r.strip().split(' ') #all fragments
                prev = r[:self.start]  #read name
                if is_contacts:
                    frags = []
                    contacts = r[self.start:]
                    for con in contacts:
                        con = con.split(',')
                        frags.extend([con[0][1:],con[1][:-1]])
                    frags = set(frags)
                else:
                    frags = set(r[self.start:])
                if frags&region==set():
                    continue
                fout.write(r_org)
    def find_region_only(self,ffout,region):
        l=0
        with open(ffout,'w') as fout:
            for r in self.fin:
                l += 1
                if l % 10000 == 0:
                    sys.stdout.write("\r%d" % l)
                    sys.stdout.flush()
                r_org = r
                r = r.strip().split(' ')
                bc = r[0]
                prev = r[:self.start]
                frags = set(r[self.start:])
                overlaps = frags&region
                if overlaps==set():
                    continue
                fout.write(bc+' '+' '.join(overlaps)+'\n')   
    def filt_sequential(self,ffout,cut=3):
        l=0
        self.reset_file()
        with open(ffout,'w') as fout:
            for r in self.fin:
                l += 1
                if l%1000000==0:
                    sys.stdout.write("\r%d" % l)
                    sys.stdout.flush()
                #if l>10000:
                #    break
                r_org = r
                r = r.strip().split(' ')
                r = r[self.start:]
                if len(r)<2:
                    continue
                ch_locs = [i.split('#') for i in r]
                chs = set([i[0] for i in ch_locs])
                if len(chs)>1:
                    fout.write(r_org)
                    continue
                locs = np.array(sorted([int(i[1]) for i in ch_locs]))
                chas = locs[1:] - locs[:-1]
                if max(chas)>=cut:
                    fout.write(r_org)