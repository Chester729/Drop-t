# Drop-t
algorithm code for Drop-t

Run the code to obtain d-LHCC from Drop-t sequencing data.
pipeline:
raw sequencing reads
use longranger to get barcoded reads
use fastp filt quality
use juicer get merged_nodups.txt
use merged_nodups to get bc_contacts
bc_contacts > bc_contacts_sorted
from bc_contacts_sorted get all bc fragments, bc_contacts


input bc_frags, bc_contacts
make_weighted_graph.py
