"""
Due to the issue of repeat region, we remove some fragments (Default 4000) appeared the most times.
We record these "top frags" for each barcode.
Keys: bc
Values: dict:
        keys:top frags, 
        values: edges of the frag
"""
from requirements.Bc_analysis import bc_frag_analysis
import sys
import argparse
from argcomplete import autocomplete
parser = argparse.ArgumentParser(description="Process some parameters.")
parser.add_argument('-f','--frag_file', required=True,help="frag file")
parser.add_argument('-c','--contact_file', required=True, help="contact file")
parser.add_argument('-w','--write', default=True, type=bool,help="genome")
parser.add_argument('-ft','--file_top', default='top_frags.txt',help="top frags fout")
parser.add_argument('-t','--top_frag_dict_out', default='top_frag_dict.npy',help="top frag dict out")
parser.add_argument('-n','--top_num',default=4000,type=int,help="top frag dict out")
autocomplete(parser)
args = parser.parse_args()
bc_frag = bc_frag_analysis(frag_file=args.frag_file,contact_file=args.contact_file)
if args.write==True:
    bc_frag.make_num_bc_frag()
    bc_frag.get_top_frags(top_num=args.top_num,write=args.write,ffout = args.file_top)
else:
    bc_frag.make_num_bc_frag()
    bc_frag.get_top_frags(file=args.file_top)
bc_frag.make_top_frag_dict(args.top_frag_dict_out)