#this is a pipeline for the examples
#learn the format of all files by running this and check the files in "examples" directory

#run XXX.py -h for usage, some parameters could be adjusted.
#make a graph
python3 make_weighted_graph_big.py -f examples/frags_10000.txt -c examples/contacts_10000.txt -o examples/test_graph.GraphML -m ./hg19_MboI.txt

#get all top fragments (usually repeats), should be removed in the following analysis
python3 make_top_frag_dict.py -f examples/frags_10000.txt -c examples/contacts_10000.txt -t examples/top_frag_dict.npy -n 400


#Directly find components, for parameter determination
python3 Find_components-delete_top.py -g examples/test_graph.GraphML -t examples/top_frag_dict.npy -f examples/frags_10000.txt -c examples/Components_total-delete_top


#parameter determination
python3 calculate_no_neighbor_prob.py -t top_frags.txt -m hg19_MboI.txt -g examples/test_graph.GraphML -c examples/Components_total-delete_top


#connection, get d-LHCCs
python3 Components_combine_final.py -m hg19_MboI.txt -g examples/test_graph.GraphML -t examples/top_frag_dict.npy -d examples/frags_10000.txt
