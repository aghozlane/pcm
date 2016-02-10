#!/bin/bash
PBS_SCRIPT=$HOME/modeller_submission.sh
tmpdir=/pasteur/scratch/amine/
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
header="""#!/bin/bash
#$ -S /bin/bash
#$ -M amine.ghozlane@pasteur.fr
##$ -m bea
#$ -q hubbioit
#$ -pe thread $CPU
#$ -l mem_total=30G
### LIBRARY

export PYTHONPATH=$HOME/soft/python-lib/lib/python2.7/site-packages/requests-2.3.0-py2.7.egg:$HOME/bin/modeller9.15/modlib/:$HOME/bin/modeller9.15/lib/x86_64-intel8/python2.5
export LD_LIBRARY_PATH=$HOME/bin/modeller9.15/lib/x86_64-intel8:$LD_LIBRARY_PATH
export prodir=$HOME/soft/procheck
export PROQ_DIR12=$HOME/soft/ProQv1.2
export NACCESS=$HOME/soft/naccess2.1.1/naccess
export STRIDE=$HOME/soft/stride/stride
export TMPDIR=$tmpdir/tmp
#export BLASTMAT=/usr/local/bioinfo/src/NCBI_Blast/current/data

module add blast/2.2.26 blast+/2.2.31 Python/2.7.8 Clustal-Omega/1.2.1
"""
# BUILD PBS scripts
for file in `ls $1/*.fasta`
do
  echo "Construction du script PBS"
  SampleFile=$(basename $file)
  SampleName="${SampleFile%.*}"
  output_dir="$dir/$SampleName/"
  output_dir2="$dir2/$SampleName/"
  modeller_summary=$(ls -1 $output_dir/modeller_summary_* 2> /dev/null | head -n 1)
  modeller_summary2=$(ls -1 $output_dir2/modeller_summary_* 2> /dev/null | head -n 1)
  if [ ! -f "$modeller_summary" ]
  then
        nb_pdb=$(ls -1 $output_dir/*.pdb 2> /dev/null |wc -l )
        mkdir -p $output_dir
        if [ "$nb_pdb" -eq 0 ]
        then
            echo "run on $output_dir"
            echo """$header
#$ -o $tmpdir/rapport/${SampleName}_ref_output.txt
#$ -e $tmpdir/rapport/${SampleName}_ref_error.txt
#$ -N "ref_${SampleName}"
python $HOME/Modeller/modeller_script.py -l model check -s prosa proq_standalone  -f $file -pi blastp -pd $4  -pr $HOME/database/cleaned_pdb/ -j $HOME/soft/psipred/ -r $output_dir  -k $HOME/soft/ProQv1.2/ -q max -n 100  -nb 3 -t $CPU &> $output_dir/log_modeller.txt""">$PBS_SCRIPT
            PBSID=`qsub $PBS_SCRIPT`
            #cat $PBS_SCRIPT
            echo "! Soumission PBS :> JOBID = $PBSID"
        else
            echo "Already done on : $output_dir"
        fi
  fi
  if [ ! -f "$modeller_summary2" ]
  then
        nb_pdb=$(ls -1 $output_dir2/*.pdb 2> /dev/null |wc -l)
        mkdir -p $output_dir2
        if [ "$nb_pdb" -eq 0 ]
        then
            echo "run on $output_dir2"
            echo """$header
#$ -o $tmpdir/rapport/${SampleName}_tneg_output.txt
#$ -e $tmpdir/rapport/${SampleName}_tneg_error.txt
#$ -N "tneg_${SampleName}"
python $HOME/Modeller/modeller_script.py -l model check -s prosa proq_standalone  -f $file -pi blastp  -pd $5  -pr $HOME/database/cleaned_pdb/ -j $HOME/soft/psipred/ -r $output_dir2  -k $HOME/soft/ProQv1.2/ -q max -n 100 -nb 3  -t $CPU &> $output_dir2/log_modeller.txt""">$PBS_SCRIPT
            PBSID=`qsub $PBS_SCRIPT`
            #cat $PBS_SCRIPT
            echo "! Soumission PBS :> JOBID = $PBSID"
       else
            echo "Already done on : $output_dir2"
       fi
  fi
done
