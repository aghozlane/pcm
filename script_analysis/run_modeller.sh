#!/bin/bash
PBS_SCRIPT=/projet/externe/inra/aghozlane/evotar/modeller_submission.sh

#Check arguments
if [ $# -ne 6  ]
then
    echo "$0 <fasta_dir> <modelling_dir_1> <modelling_dir_2> <database_1> <database_2> <nb_cpu>"
    exit
fi
#Outdir
dir="$(pwd)/$2"
dir2="$(pwd)/$3"
CPU=$6
header="""#!/bin/bash
#$ -S /bin/bash
#$ -M amine.ghozlane@jouy.inra.fr
#$ -m bea
#$ -q short.q
#$ -pe thread $CPU
### LIBRARY
export PYTHONPATH=$HOME/soft/python-lib/lib/python2.7/site-packages/requests-2.3.0-py2.7.egg:$HOME/bin/modeller9.14/modlib/:$HOME/bin/modeller9.14/lib/x86_64-intel8/python2.5
export LD_LIBRARY_PATH=$HOME/bin/modeller9.14/lib/x86_64-intel8:$LD_LIBRARY_PATH
export prodir=$HOME/soft/procheck
export PROQ_DIR12=/home/externe/inra/aghozlane/soft/ProQv1.2/
export NACCESS=/home/externe/inra/aghozlane/soft/naccess2.1.1/naccess
export STRIDE=/home/externe/inra/aghozlane/soft/stride/stride
export TMPDIR=/projet/externe/inra/aghozlane/evotar/tmp
"""
# BUILD PBS scripts
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
        mkdir -p $output_dir
        if [ "$nb_pdb" -eq 0 ]
        then
            echo "run on $output_dir"
            echo """$header
#$ -o /projet/externe/inra/aghozlane/evotar/tmp/${SampleName}_ref_output.txt
#$ -e /projet/externe/inra/aghozlane/evotar/tmp/${SampleName}_ref_error.txt
#$ -N "ref_${SampleName}"
python /projet/externe/inra/aghozlane/evotar/script/Modeller/modeller_script.py -l model check -s prosa proq_standalone  -f $(pwd)/$file -pi blastp -pd $(pwd)/$4  -pr /projet/externe/inra/aghozlane/evotar/input/database/cleaned_pdb/ -j $HOME/soft/psipred/ -r $output_dir  -k  $HOME/soft/ProQv1.2/ -q max -n 100  -nb 3 -t $CPU -g /home/externe/inra/aghozlane/soft/clustalo/ -pp /usr/local/genome2/bin/ &> $output_dir/log_modeller.txt""">$PBS_SCRIPT
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
#$ -o /projet/externe/inra/aghozlane/evotar/tmp/${SampleName}_tneg_output.txt
#$ -e /projet/externe/inra/aghozlane/evotar/tmp/${SampleName}_tneg_error.txt
#$ -N "tneg_${SampleName}"
python /projet/externe/inra/aghozlane/evotar/script/Modeller/modeller_script.py -l model check -s prosa proq_standalone  -f $(pwd)/$file -pi blastp  -pd $(pwd)/$5  -pr /projet/externe/inra/aghozlane/evotar/input/database/cleaned_pdb/ -j $HOME/soft/psipred/ -r $output_dir2  -k  $HOME/soft/ProQv1.2/ -q max -n 100 -nb 3 -t $CPU -g /home/externe/inra/aghozlane/soft/clustalo/ -pp /usr/local/genome2/bin/ &> $output_dir2/log_modeller.txt""">$PBS_SCRIPT
            PBSID=`qsub $PBS_SCRIPT`
            #cat $PBS_SCRIPT
            echo "! Soumission PBS :> JOBID = $PBSID"
       else
            echo "Already done on : $output_dir2"
       fi
  fi
done
