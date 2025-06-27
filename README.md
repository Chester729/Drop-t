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

Required files

bc_frags: all restriction fragments in each droplet

format: (See file examples/frags_10000.txt)

bc1  frag11  frag12  frag13  ...

bc2  frag21  frag22  frag23  ...

...

bc_contacts: all contacts between fragments in each droplet

format: (See file examples/contacts_10000.txt)

bc1  (frag111,frag112):weight  (frag121,frag122):weight  ...

bc2  (frag111,frag212):weight  (frag221,frag222):weight  ...

...

MboI.txt: genomic coordinates of all MboI (or other enzyme) restriction site 
This can be obtained from juicer.

codes:
make_weighted_graph.py: from bc_frags and bc_contacts create a weighted graph

