#!/bin/bash

function error {
    # echo red
    echo -e  "\e[31m* $1\e[0m" >&2
}



#Check arguments
if [ $# -ne 1  ]
then
     echo "$0 <result_dir>"
     exit
fi

# Check if ref and tneg dir exist
if [ ! -d "$1/ref" ] || [ ! -d "$1/tneg" ]
then
    error "$1/ref or $1/tneg doesn't exist"
fi

datasets=($1/ref $1/tneg)
project_name=$(basename $1)
results=($1/${project_name}_by_ref.csv $1/${project_name}_by_tneg.csv)

###BLAA2
#datasets=(/export/mgps/home/aghozlane/data_evotar/three_templates/blaa2/blaa2_candidates_by_ref /export/mgps/home/aghozlane/data_evotar/three_templates/blaa2/blaa2_candidates_by_tneg)
#results=(blaa2_candidates_by_ref.csv blaa2_candidates_by_tneg.csv)

#one shot
#datasets=(/datatank/transit/data_eruppe/blad/blad_ref_pcm/blad_ref/SPEA)
#results=(resume_SPEA.csv)

modeller="$(dirname "${BASH_SOURCE[0]}")/../modeller_script.py"
for i in `seq 0 1`
#for i in `seq 0 0`
do
    echo "Ecriture du fichier : ${results[$i]}"
    echo -e "Model\tmolpdf\tdope_score\tnormalized_dope\tGA341_score\tzscore\tmaxsub\tlgscore"  > ${results[$i]}
    for rep in $(ls -d ${datasets[$i]}/*/)
    do 
        #find "$rep" -name "modeller_summary_*.csv" |head -1
        summary=$(ls -1 $rep/modeller_summary_*.csv  2>/dev/null |head -1)
        #summary=$(find "$rep" -name "modeller_summary_*.csv" |head -1)
        prosa=$(ls -1 $rep/result_prosa_*.txt  2>/dev/null |head -1)
        #prosa=$(find "$rep" -name "/result_prosa_*.txt"|head -1)
        proq=$(ls -1 $rep/result_proq_*.txt  2>/dev/null |head -1 )
        #proq=$(find "$rep" -name "/result_proq_*.txt"|head -1)
        horiz=$(ls -1 $rep/*.horiz  2>/dev/null |head -1 )
        #horiz=$(find "$rep" -name "/*.horiz"|head -1)
        if [ -f "$summary" ]
        then
            if [ ! -f "$prosa" ] || [ ! -f "$proq" ]
            then
                echo "PROSA file is missing for $i, start to re-run modeller_script" >&2
                echo "$summary">&2
                $modeller -s proq_standalone prosa -l check -sm $summary -d $horiz -k $HOME/soft/ProQv1.2/   1>&2  # 2> $rep/log_check.txt  
                prosa=$(ls -1 $rep/result_prosa_*.txt  2>/dev/null |head -1)
                proq=$(ls -1 $rep/result_proq_*.txt  2>/dev/null |head -1)
                horiz=$(ls -1 $rep/*.horiz  2>/dev/null |head -1 )
                echo "done..." >&2
            fi
            if [  -f  "$prosa"  ]
            then 
               best_model=$(tail -n +2 $summary |head -1 |cut -s -f1|sed -e  's/\r//g'|sed 's/\.pdb//g')
               dope=$(tail -n +2 $summary |head -1 |cut -s -f3|sed -e  's/\r//g')
               molpdf=$(tail -n +2 $summary |head -1 |cut -s -f2|sed -e  's/\r//g')
               normalized_dope=$(tail -n +2 $summary |head -1 |cut -s -f4|sed -e  's/\r//g')
               GA341_score=$(tail -n +2 $summary |head -1 |cut -s -f5|sed -e  's/\r//g')
               zscore=$(tail -n +2 $prosa |head -1 |awk '{print $2}'|sed -e  's/\r//g'  )
               maxsub=$(tail -n +2 $proq |head -1 |awk '{print $2}'| sed -e  's/\r//g' )
               lgscore=$(tail -n +2 $proq |head -1 |awk '{print $3}'| sed -e  's/\r//g')
               echo -e "$best_model\t$molpdf\t$dope\t$normalized_dope\t$GA341_score\t$zscore\t$maxsub\t$lgscore"
            fi
        else
            echo "Modeller_summary file is missing for $rep" >&2
        fi
    done   >> ${results[$i]}
done
