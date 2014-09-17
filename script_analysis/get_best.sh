#!/bin/bash

if [[ $# -ne 2 ]]
then 
   echo "Usage : $0 ./path/to/interest/directory/ output/directory/"
   exit 1
fi

for i in `ls -d $1/*/` 
do
	pdb_model=$(head -n 2 $i/modeller_summary_*.csv |tail -n 1 |cut -f1 -s)
	echo $pdb_model
	cp $i/$pdb_model $2/
done
