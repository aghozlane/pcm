#!/bin/bash

function check_dir {
    # Check resultDir
    if [ ! -d $resultDir ]
    then
        error "$1 does not exist !"
    fi
}


#Check arguments
if [ $# -ne 3  ]
then
    echo "$0 <candidate_modelling_dir> <reference_dir> <output_dir>"
    exit 1
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
        echo "$best_model"
        cp $rep/$best_model $candidate_dir/
        #fuck=$(basename $rep)
        #fuck_suite=${fuck:3:${#fuck}}
        #cp $rep/$best_model $candidate_dir/${fuck}${fuck_suite}.pdb
    fi
done


# Select reference ?

# Get reference ?

# run alignment
script=/datatank/transit/data_eruppe/pcm
$script/PDBRMSD/PDBRMSD.py -q $candidate_dir/ -t $2/ -s TMalign mammoth  -p $script/soft/tmscore/ $script/soft/mammoth_compila/ -r $3/ -b 1
