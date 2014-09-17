#!/bin/bash

function check_dir {
    # Check resultDir
    if [ ! -d $resultDir ]
    then
        error "$1 does not exist !"
    fi
}


#Check arguments
if [ $# -ne 4 ]
then
    echo "$0 <candidate_modelling_dir> <reference_dir_noNA> <output_dir> <reference_position_noNA>"
fi

# Check resultDir
if [ ! -d $3 ]
then
        mkdir $3
        check_dir $3
fi

# Get best models
candidate_dir=$3/query/
if [ ! -d $candidate_dir ]
then
    mkdir $candidate_dir
    check_dir $3
fi

# Copy best model in candidate_dir
for rep in `ls -d $1/*/`
do
    summary=$(ls $rep/modeller_summary_*.csv 2>/dev/null)
    if [ -f "$summary" ]
    then
        best_model=$(tail -n +2 $summary |head -1 |cut -s -f1|sed -e  's/\r//g')
        fuck=$(basename $rep)
        fuck_suite=${fuck:3:${#fuck}}
        cp $rep/$best_model $candidate_dir/${fuck}${fuck_suite}.pdb
    fi
done


# run alignment
$HOME/PDBRMSD/PDBRMSD.py -q $candidate_dir/ -t $2/ -s TMalign  -p $HOME/soft/tmscore/ -r $3/ -b 10 -i $4
