# Drop-t
algorithm code for Drop-t


Before running these codes, you need to get all fragments and contacts in each droplet(barcode). 

Detailed pipeline can be find in Methods and Supplementary methods.

From raw sequencing reads,

use longranger to get barcoded reads

use fastp filter low quality reads and adapters

use juicer get merged_nodups.txt

from merged_nodups.txt get the bc_frags and bc_contacts


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

To find the order of running these codes, check "example_codes.sh". 

To find what each code does, check example_codes.sh and the .py files.
