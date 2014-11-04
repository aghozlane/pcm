#!/bin/bash


#Check arguments
if [ $# -ne 6  ]
then
    echo "$0 <fasta_dir> <modelling_dir_1> <modelling_dir_2> <database_1> <database_2> <nb_cpu>"
    exit
fi

#Outdir
dir="$2"
dir2="$3"
CPU=$6

for file in `ls $1/*.fasta`
do
  echo "Contruction du script PBS"
  SampleFile=$(basename $file)
  SampleName="${SampleFile%.*}"
  output_dir="$dir/$SampleName/"
  output_dir2="$dir2/$SampleName/"
  modeller_summary=$(ls -1 $output_dir/modeller_summary_* 2> /dev/null | head -n 1)
  modeller_summary2=$(ls -1 $output_dir2/modeller_summary_* 2> /dev/null | head -n 1)
  if [ ! -f "$modeller_summary" ]
  then
        nb_pdb=$(ls -1 $output_dir/*.pdb 2> /dev/null |wc -l )
        #mkdir -p $output_dir
        if [ "$nb_pdb" -eq 0 ]
        then
            echo "run on $output_dir"
            python $HOME/Modeller/modeller_script.py -l model check -s prosa proq_standalone  -f $file -pi blastp -pd $4  -pr $HOME/result/metahit_etienne/cleaned_pdb/ -j $HOME/soft/psipred/ -r $output_dir  -k $HOME/soft/procheck/ $HOME/soft/ProQv1.2/ -q max -n 100  -nb 3 -t $CPU &> $output_dir/log_modeller.txt
        else
            echo "Already done on : $output_dir"
        fi
  fi
  if [ ! -f "$modeller_summary2" ]
  then
        nb_pdb=$(ls -1 $output_dir2/*.pdb 2> /dev/null |wc -l)
        #mkdir -p $output_dir2
        if [ "$nb_pdb" -eq 0 ]
        then
            echo "run on $output_dir2"
            python $HOME/save/Modeller/modeller_script.py -l model check -s prosa proq_standalone  -f $file -pi blastp  -pd $5 -pr $HOME/result/metahit_etienne/cleaned_pdb/ -j $HOME/soft/psipred/ -r $output_dir2  -k $HOME/soft/procheck/ $HOME/soft/ProQv1.2/ -q max -n 100 -nb 3  -t $CPU &> $output_dir2/log_modeller.txt
       else
            echo "Already done on : $output_dir2"
       fi
  fi
done

