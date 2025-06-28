import gzip
import numpy as np
import sys
from .make_MboI import make_MboI
class bc_frag_analysis:
    def __init__(self, frag_file=None, contact_file=None,genome='hg19',get_MboI=False,add_chr=True):
        self.genome = genome
        self.frag_file=frag_file
        self.contact_file=contact_file
        self.add_chr=add_chr
        # Initialize MboI
        #print(self.genome,self.add_chr,line1[-1][:2].lower())
        if get_MboI:
            self.chrindex,self.MboI = make_MboI(self.genome, self.add_chr)
    def set_frag_file(self,frag_file):
        self.frag_file = frag_file
    def set_contact_file(self,contact_file):
        self.contact_file = contact_file
    def reset_frag_file(self):
        if self.frag_file.endswith('.gz'):
            self.fin_frag = gzip.open(self.frag_file, 'rt')
        else:
            self.fin_frag = open(self.frag_file, 'r')
        
    def reset_contact_file(self):
        if self.contact_file.endswith('.gz'):
            self.fin_contact = gzip.open(self.contact_file, 'rt')
        else:
            self.fin_contact = open(self.contact_file, 'r')
    def barcode_fragnum(self):
        self.reset_frag_file()
        size=[]
        l = 0
        for r in self.fin_frag:
            l+=1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            size.append(len(r[1:]))
            #if l>10000:
            #    break
        self.barcode_size = np.array(size)
    def barcode_contactnum(self):
        self.reset_contact_file()
        size=[]
        l = 0
        for r in self.fin_contact:
            l+=1
            if l % 1000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            size.append(len(r[1:]))
            #if l>10000:
            #    break
        self.barcode_contacts = np.array(size)
    def get_add_chr(self):
        self.reset_frag_file()
        line1 = self.fin.readline().strip().split(' ')
        self.add_chr = line1[-1][:2].lower() == 'ch'
    def set_add_chr(add_chr):
        self.add_chr=add_chr
    def make_num_bc_frag(self):   #每个fragment出现次数
        num_bc_frag = {}
        l=0
        self.reset_frag_file()
        for r in self.fin_frag:
            l += 1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            bc = r[0]
            for f in r[1:]:
                if f not in num_bc_frag.keys():
                    num_bc_frag[f]=0
        print('\n',len(num_bc_frag))
        l=0
        self.reset_frag_file()
        for r in self.fin_frag:
            l += 1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            bc = r[0]
            for f in r[1:]:
                try:
                    num_bc_frag[f]+=1
                except KeyError:
                    pass
        self.num_bc_frag = num_bc_frag
    def split_frags_by_chrs(self,query_chr,ffout,ffout_other=None):
        query_chr = set(query_chr)
        self.reset_frag_file()
        fout = open(ffout,'w')
        if ffout_other!=None:
            fout_other=open(ffout_other,'w')
        l=0
        for r in self.fin_frag:
            l += 1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            frags = r[1:]
            bc = r[0]
            query_frags = [i for i in frags if i.split('#')[0] in query_chr]
            fout.write(bc+' '+' '.join(query_frags)+'\n')
            if ffout_other!=None:
                other_frags = [i for i in frags if i.split('#')[0] not in query_chr]
                fout_other.write(bc+' '+' '.join(other_frags)+'\n')
    def split_contacts_by_chrs(self,query_chr,ffout,ffout_other=None,ffout_hybrid=None):
        query_chr = set(query_chr)
        self.reset_contact_file()
        fout = open(ffout,'w')
        if ffout_other!=None:
            fout_other=open(ffout_other,'w')
        if ffout_hybrid!=None:
            fout_hybrid = open(ffout_hybrid,'w')
        l=0
        for r in self.fin_contact:
            l += 1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            bc = r[0]
            contacts = r[1:]
            query_contacts = []
            other_contacts = []
            hybrid_contacts = []
            for i in contacts:
                i_split = i.split(',')
                frag1 = i_split[0][1:]
                frag2 = i_split[1][:-1]
                ch1 = frag1.split('#')[0]
                ch2 = frag2.split('#')[0]
                query_ch_num = (ch1 in query_chr) + (ch2 in query_chr)
                if query_ch_num==2:
                    query_contacts.append(i)
                elif query_ch_num==1:
                    hybrid_contacts.append(i)
                else:
                    other_contacts.append(i)
            fout.write(bc+' '+' '.join(query_contacts)+'\n')
            if ffout_other!=None:
                fout_other.write(bc+' '+' '.join(other_contacts)+'\n')               
            if ffout_hybrid!=None:
                fout_hybrid.write(bc+' '+' '.join(hybrid_contacts)+'\n')   
    def write_top_frags(self,ffout):
        with open(ffout,'w') as fout: 
            for f in self.top_frags:
                fout.write(f+' ')
    def get_top_frags(self,top_num=4000,write=False,ffout=None,file=None):
        if file==None:
            self.top_frags = [i[0] for i in sorted(self.num_bc_frag.items(),key=lambda x:-x[1])[:top_num]]
        else:
            with open(file,'r') as fin:
                self.top_frags = set(fin.readline().split(' '))
        if write:
            self.write_top_frags(ffout)
    def make_top_frag_dict(self,save_file=None):
        top_frag_dict = {}
        self.reset_frag_file()
        self.reset_contact_file()
        l=0
        while True:
            l += 1
            if l % 10000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            frags = self.fin_frag.readline().split(' ')
            contacts = self.fin_contact.readline().split(' ')
            frags[-1] = frags[-1][:-1]
            contacts[-1] = contacts[-1][:-1]
            if frags[0]=='' or contacts[0]=='':
                break
            bc = frags[0]
            bc2 = contacts[0]
            if bc!=bc2:
                print(bc,bc2)
                print('different barcode!\n')
                return 0
            top_frag_dict[bc] = {}
            frags = frags[1:]
            contacts = contacts[1:]
            overlaps = set(frags)&set(self.top_frags)
            for i in overlaps:
                top_frag_dict[bc][i] = []
            for i in contacts:
                pairs = i.split(':')
                if i=='':
                    continue
                weight = int(pairs[1])
                pairs = pairs[0]
                pairs = pairs.split(',')
                frag1 = pairs[0][1:]
                frag2 = pairs[1][:-1]
                if frag1 in overlaps:
                    top_frag_dict[bc][frag1].append((frag1,frag2))
                elif frag2 in overlaps:
                    top_frag_dict[bc][frag2].append((frag1,frag2))
        self.top_frag_dict = top_frag_dict
        if save_file!=None:
            np.save(save_file,self.top_frag_dict)
    ###实际切割大小
    def get_sequential_length(self,interval=3):
        self.reset_frag_file()
        l=0
        core_lens = []
        for r in self.fin_frag:
            l+=1
            if l % 1000 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            r = r.strip().split(' ')
            r = r[1:]
            ch_loc = [j.split('#') for j in r]
                #print(ch_loc)
            ch_loc = [(j[0],int(j[1])) for j in ch_loc]
            ch_loc = sorted(ch_loc)
            now_core = [ch_loc[0]]
            for j in range(len(ch_loc)-1):
                if ch_loc[j][0]==ch_loc[j+1][0]:   #same chr
                    now_gap = ch_loc[j+1][1]-ch_loc[j][1]
                    if now_gap<interval:
                        now_core.append(ch_loc[j+1])
                    else:
                        try:
                            s=0
                            for i in now_core:
                                s+=self.MboI[i[0]][i[1]+1]-self.MboI[i[0]][i[1]]
                        except IndexError:
                            break
                        core_lens.append(s)
                        now_core = [ch_loc[j+1]]
                else:
                    try:
                        s=0
                        for i in now_core:
                            s+=self.MboI[i[0]][i[1]+1]-self.MboI[i[0]][i[1]]
                    except IndexError:
                        break
                    now_core = [ch_loc[j+1]]
                    core_lens.append(s)
                #if l>0:
                #    break
        self.core_lens = core_lens
