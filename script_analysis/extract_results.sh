#!/bin/bash

#check arguments
if [ $# -ne "4" ]
then
    echo "Usage : $0 <reference_dir> <negative_dir> <reference_result> <negative_result>"
    exit
fi



datasets=("$1" "$2")
results=("$3" "$4")


script="$HOME/Modeller/"
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
        if [ ! -f "$prosa" ] 
        then
            echo "PROSA file is missing for $i, start to re-run modeller_script" >&2
            #echo "$script/Modeller/modeller_script.py -k $script/soft/procheck/ -l check -sm $summary -d $horiz"
            echo "$summary">&2
            $script/Modeller/modeller_script.py -s prosa -l check -sm $summary -d $horiz 1>&2  # 2> $rep/log_check.txt  
            prosa=$(ls -1 $rep/result_prosa_*.txt  2>/dev/null |head -1)
            proq=$(ls -1 $rep/result_proq_*.txt  2>/dev/null |head -1)
            horiz=$(ls -1 $rep/*.horiz  2>/dev/null |head -1 )
            echo "done..." >&2
        fi
        if [ ! -f "$proq" ]
            echo "proq file is missing for $i, start to re-run modeller_script" >&2
            #echo "$script/Modeller/modeller_script.py -k $script/soft/procheck/ -l check -sm $summary -d $horiz"
            echo "$summary">&2
            $script/Modeller/modeller_script.py -s proq_standalone -l check -sm $summary -d $horiz 1>&2  # 2> $rep/log_check.txt  
            prosa=$(ls -1 $rep/result_prosa_*.txt  2>/dev/null |head -1)
            proq=$(ls -1 $rep/result_proq_*.txt  2>/dev/null |head -1)
            horiz=$(ls -1 $rep/*.horiz  2>/dev/null |head -1 )
            echo "done..." >&2
        fi
        if [  -f  "$prosa"  ]
        then 
            best_model=$(tail -n +2 $summary |head -1 |cut -s -f1|sed -e  's/\r//g')
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
