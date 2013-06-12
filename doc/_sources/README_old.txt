Model_laccase script

================
=1. QUICK START=
================
==================
1.1. Requirements
==================
- On IENA : 
module add modeller/9.11

- Otherwise add the following expression in your bashrc :
export PYTHONPATH=/usr/lib/modeller9.11/modlib:/usr/lib/modeller9.11/lib/x86_64-intel8/python2.5:$PYTHONPATH
export LD_LIBRARY_PATH=/usr/lib/modeller9.11/lib/x86_64-intel8:$LD_LIBRARY_PATH

And :
source .bashrc

==================
1.2. Run
==================
To get the description of script options :

modeller_script.py -h

- First case: one alignment (template + model) + pdb (template)

modeller_script.py  -i alignment.pir  -p structure.pdb

- Second case: one multifasta (template + model) + pdb (template) + alignment with clustalo

modeller_script.py  -f sequences.fasta  -p structure.pdb  -a clustalo -g path/to/clustalo/

- Third case: modelize with several template and set 3 hetero atom from only one of these structures:

modeller_script.py   -f sequences.fasta   -p structure1.pdb structure2.pdb -ht 3 -hm Model_name Structure1_name
