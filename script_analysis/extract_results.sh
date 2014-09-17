#!/bin/bash

###BLAA
#datasets=(/datatank/transit/data_aghozlane/three_templates/blaa_forsberg_complement/blaa_forsberg_complement_ref /datatank/transit/data_aghozlane/three_templates/blaa_forsberg_complement/blaa_forsberg_complement_tneg)
#results=(blaa_forsberg_complement_by_ref.csv blaa_forsberg_complement_by_tneg.csv)

###BLAA2
datasets=(/export/mgps/home/aghozlane/data_evotar/three_templates/blaa2/blaa2_candidates_by_ref /export/mgps/home/aghozlane/data_evotar/three_templates/blaa2/blaa2_candidates_by_tneg)
results=(blaa2_candidates_by_ref.csv blaa2_candidates_by_tneg.csv)

#one shot
#datasets=(/datatank/transit/data_eruppe/blad/blad_ref_pcm/blad_ref/SPEA)
#results=(resume_SPEA.csv)

script="/datatank/transit/data_eruppe/pcm"
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
            #echo "$script/Modeller/modeller_script.py -k $script/soft/procheck/ -l check -sm $summary -d $horiz"
            echo "$summary">&2
            $script/Modeller/modeller_script.py -s proq prosa verify3D -l check -sm $summary -d $horiz 1>&2  # 2> $rep/log_check.txt  
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
            
            #fuck=$(basename $rep)
            #fuck_suite=${fuck:3:${#fuck}}

            echo -e "$best_model\t$molpdf\t$dope\t$normalized_dope\t$GA341_score\t$zscore\t$maxsub\t$lgscore"
            #echo -e "${fuck}${fuck_suite}\t$molpdf\t$dope\t$normalized_dope\t$GA341_score\t$zscore\t$maxsub\t$lgscore"
        fi
    else
        echo "Modeller_summary file is missing for $rep" >&2
    fi
done   >> ${results[$i]}
done
