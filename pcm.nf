#!/usr/bin/env nextflow

/*workflow.onComplete {
    def subject = 'Annot virus'
    def recipient = 'amine.ghozlane@pasteur.fr'

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}*/

workflow.onComplete = {
    // any workflow property can be used here
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}


workflow.onError = {
    println "Oops .. something went wrong"
}

// General parameters
params.in = "${baseDir}/test/"
params.out = "${baseDir}/test/res"
params.modelling = "${params.out}/modelling"
params.cpu = 6
//params.database = "${baseDir}/database/"
params.database = "/usr/local/bin/database/"
params.universal_model = "${params.database}/universal_model.csv"
params.cleaned_pdb = "${params.database}/cleaned_pdb/"
//params.proq = "$HOME/soft/ProQv1.2/"
params.proq = "/usr/local/bin/"
//params.psipred = "$HOME/soft/psipred/"
params.psipred = "/usr/local/bin/"
params.internet_access = "no"
params.modelling_quality = "fast"
params.num_model = 6
params.num_template = 3
params.hfinder_evalue = 1E-5
params.number_bootstrap = 10
//params.tmalign_dir = "${baseDir}/soft/"
params.tmalign_dir = "/usr/local/bin/"
//params.mammoth_dir = "${baseDir}/soft/"
params.mammoth_dir = "/usr/local/bin/"

myDir = file(params.out)
myDir.mkdirs()
modDir = file(params.modelling)
modDir.mkdirs()

params.help=false

def usage() {
    println("pcm.nf --in <fasta_file> --out <output_dir> --cpus <nb_cpus> -w <temp_work_dir>")
    println("--in Multifasta file containing protein sequence (default ${params.in}).")
    println("--out Output directory (default ${params.out}). ")
    println("--cpu Number of cpus for process (default ${params.cpu})")
}


if(params.help){
    usage()
    exit(1)
}

multifastaChannel = Channel
                    .fromPath("${params.in}")
                    .ifEmpty { exit 1, "Missing parameter: ${params.in}" }

familyChannel = Channel
                 .from(
                    ["aac2", "142", "236", "${params.database}/aac2/aac2_ref.faa", "${params.database}/aac2/aac2.hmm"],
                    ["aac3_1", "122", "203","${params.database}/aac3_1/aac3_1_ref.faa","${params.database}/aac3_1/aac3_1.hmm"],
                    ["aac3_2", "205", "341","${params.database}/aac3_2/aac3_2_ref.faa","${params.database}/aac3_2/aac3_2.hmm"],
                    ["aac6", "123", "205","${params.database}/aac6/aac6_ref.faa","${params.database}/aac6/aac6.hmm"],
                    ["ant", "194", "323","${params.database}/ant/ant_ref.faa","${params.database}/ant/ant.hmm"],
                    ["aph", "206", "344","${params.database}/aph/aph_ref.faa","${params.database}/aph/aph.hmm"],
                    ["arnm", "190", "316","${params.database}/arnm/arnm_ref.faa","${params.database}/arnm/arnm.hmm"],
                    ["blaa", "219", "365","${params.database}/blaa/blaa_ref.faa","${params.database}/blaa/blaa.hmm"],
                    ["blab1", "191", "318","${params.database}/blab1/blab1_ref.faa","${params.database}/blab1/blab1.hmm"],
                    ["blab3", "215", "359","${params.database}/blab3/blab3_ref.faa","${params.database}/blab3/blab3.hmm"],
                    ["blac", "290", "483","${params.database}/blac/blac_ref.faa","${params.database}/blac/blac.hmm"],
                    ["blad", "203", "338","${params.database}/blad/blad_ref.faa","${params.database}/blad/blad.hmm"],
                    ["dfra", "114", "236","${params.database}/dfra/dfra_ref.faa","${params.database}/dfra/dfra.hmm"],
                    ["erm", "203", "338","${params.database}/erm/erm_ref.faa","${params.database}/erm/erm.hmm"],
                    ["fos", "105", "175","${params.database}/fos/fos_ref.faa","${params.database}/fos/fos.hmm"],
                    ["ldt", "257", "429","${params.database}/ldt/ldt_ref.faa","${params.database}/ldt/ldt.hmm"],
                    ["mcr", "172", "578","${params.database}/mcr/mcr_ref.faa","${params.database}/mcr/mcr.hmm"],
                    ["qnr", "164", "274","${params.database}/qnr/qnr_ref.faa","${params.database}/qnr/qnr.hmm"],
                    ["sul", "209", "349","${params.database}/sul/sul_ref.faa","${params.database}/sul/sul.hmm"],
                    ["tetM", "475", "791","${params.database}/tetM/tetM_ref.faa","${params.database}/tetM/tetM.hmm"],
                    ["tetX", "287", "478","${params.database}/tetX/tetX_ref.faa","${params.database}/tetX/tetX.hmm"],
                    ["van", "260", "433","${params.database}/van/van_ref.faa","${params.database}/van/van.hmm"]
                    )

// index
process index_query {
    input:
    file(fasta) from multifastaChannel

    output:
    set file(fasta), file("*") into multifastaindexedChannel

    shell:
    """
    makeblastdb -in !{fasta} -dbtype prot
    """
}

process search_distant_homologuous {
    tag "${fam[0]}"
    publishDir "$myDir/candidates/", mode: 'copy'
    cpus params.cpu

    input:
    set file(fasta), file(index) from multifastaindexedChannel
    each fam from familyChannel

    output:
    set val("${fam[0]}"), file("*_candidates/*.fasta") optional true into candidatesChannel
    file("*_candidates/*.tsv") optional true into summaryChannel

    shell:
    """
    mkdir !{fam[0]}_candidates/
    hfinder.py -q !{fam[3]} -qm !{fam[4]} -d ${fasta} -s blastp hmmsearch ssearch -e !{params.hfinder_evalue} -lmin !{fam[1]} -lmax !{fam[2]} -b extract cumulative check -r !{fam[0]}_candidates/ -n 1 -t !{params.cpu}
    if [ -f "!{fam[0]}_candidates/all_protein_homology.fasta" ]
    then
        mv !{fam[0]}_candidates/all_protein_homology.fasta !{fam[0]}_candidates/!{fam[0]}_candidates.fasta
        mv !{fam[0]}_candidates/all_hit_length.tsv !{fam[0]}_candidates/!{fam[0]}_candidates_hit_length.tsv
        mv !{fam[0]}_candidates/all_hit_properties.tsv !{fam[0]}_candidates/!{fam[0]}_candidates_hit_properties.tsv
    fi
    """
}

process fastaExtract {
    tag "${fasta}"

    input:
    set fam, file(fasta) from candidatesChannel

    output:
    set fam, file("splitted/*.fasta") into fastaChannel mode flatten

    shell:
    """
    extract_sequence.py !{fasta} splitted/
    """
}



//if  internet access on nodes or not
if ( params.internet_access == "no") {
    process homology_modelling {
        tag "${fasta.baseName}:${fam}"
        publishDir "$myDir/modelling/${fam}_candidates/", mode: 'copy'
        cpus params.cpu
        label 'modelling'
        validExitStatus 0,3

        input:
        set fam, file(fasta) from fastaChannel

        output:
        set fam, file(fasta), file("*.horiz"), file("best_model_ref/*.pdb"), file("ref/*/result_proq_*"), file("ref/*/modeller_summary_*.csv"), file("best_model_tneg/*.pdb"), file("tneg/*/result_proq_*"), file("tneg/*/modeller_summary_*.csv") optional true into modelChannel
        file("*/*/*.svg") into imgChannel
        file("*/*/*.pdb") into allPDBChannel


        shell:
        """
        !{params.psipred}/runpsipred !{fasta} !{params.cpu}
        mkdir -p ref/!{fasta.baseName}/ best_model_ref/ tneg/!{fasta.baseName}/ best_model_tneg/
        # ref
        modeller_script_singularity.py -l model check -s proq_standalone -f !{fasta} -r ref/!{fasta.baseName}/ -pd !{params.database}/!{fam}/!{fam}_ref_pdb.faa -pr !{params.cleaned_pdb} -k !{params.proq} -d !{fasta.baseName}.horiz -q !{params.modelling_quality} -n !{params.num_model}  -nb !{params.num_template} -t !{params.cpu}  -pi blastp
        # tneg
        modeller_script_singularity.py -l model check -s proq_standalone -f !{fasta} -r tneg/!{fasta.baseName}/ -pd !{params.database}/!{fam}/!{fam}_tneg_pdb.faa -pr !{params.cleaned_pdb} -k !{params.proq} -d !{fasta.baseName}.horiz -q !{params.modelling_quality} -n !{params.num_model}  -nb !{params.num_template} -t !{params.cpu}  -pi blastp
        # Get the best model for ref
        summary_ref=\$(ls -1 ref/!{fasta.baseName}//modeller_summary_*.csv  2>/dev/null |head -1)
        # Check ref file
        if [ -f "\$summary_ref" ]
        then
            cp ref/!{fasta.baseName}//\$(sed -n 2p \$summary_ref |cut -f 1) best_model_ref/
        else
            echo "\$summary_ref file is missing"
        fi
        # Get the best model for tneg
        summary_tneg=\$(ls -1 tneg/!{fasta.baseName}//modeller_summary_*.csv  2>/dev/null |head -1)
        # Check tneg file
        if [ -f "\$summary_tneg" ]
        then
            cp tneg/!{fasta.baseName}//\$(sed -n 2p \$summary_tneg |cut -f 1) best_model_tneg/
        else
            echo "\$summary_tneg file is missing"
        fi
        """
    }

    // modelChannel.subscribe{ println it }

    process prosa_check {
         tag "${fasta.baseName}:${fam}"
         publishDir "$myDir/modelling/${fam}_candidates/", mode: 'copy'
         label 'modelling'
         errorStrategy 'retry'

         input:
         set val(fam), file(fasta), horiz, best_pdb_ref, proq_ref, summary_ref, best_pdb_tneg, proq_tneg, summary_tneg from modelChannel

         output:
         set val(fam), file(fasta), best_pdb_ref, proq_ref, file("ref/*/result_prosa_*"), summary_ref, best_pdb_tneg, proq_tneg, file("tneg/*/result_prosa_*"), summary_tneg into extractChannel
         set val(fam), file(fasta), best_pdb_ref, best_pdb_tneg into structuralAlignmentChannel

         shell:
         """
         mkdir -p ref/!{fasta.baseName}/ tneg/!{fasta.baseName}/
         modeller_script_singularity.py -s prosa -l check -sm !{summary_ref} -d !{horiz}
         modeller_script_singularity.py -s prosa -l check -sm !{summary_tneg} -d !{horiz}
         cp \$(dirname !{summary_ref})/result_prosa_* ref/!{fasta.baseName}/
         cp \$(dirname !{summary_tneg})/result_prosa_* tneg/!{fasta.baseName}/
         """
    }
}

// else{
//     process compute_pcm_full {
//         tag "${fasta.baseName}:${fam}"
//         publishDir "$myDir/modelling/", mode: 'copy'
//         cpus params.cpu
//         label 'modelling'

//         input:
//         set fam, fasta from fastaChannel

//         output:
//         set file("ref/*/*.pdb"), file("ref/*/result_*"), file("ref/*/modeller_summary_*.csv"), file("ref/*/*.horiz") into refChannel
//         set file("ref/*/*.png"), file("ref/*/*.svg") into imgChannel

//         shell:
//         """
//         mkdir -p ref/!{fasta.baseName}/
//         #cd ref/!{fasta.baseName}/
//         modeller_script_singularity.py -l model check -s proq_standalone prosa -f !{fasta} -r ref/!{fasta.baseName}/ -pd !{params.database}/!{fam}/!{fam}_ref_pdb.faa -pr !{params.cleaned_pdb} -k !{params.proq} -j !{params.psipred} -q max -n 100  -nb 3 -t !{params.cpu}  -pi blastp
//         """
//     }
// }

process extract_result {
     tag "${fasta.baseName}:${fam}"

     input:
     set val(fam), file(fasta), best_pdb_ref, proq_ref, prosa_ref, summary_ref, best_pdb_tneg, proq_tneg, prosa_tneg, summary_tneg from extractChannel

     output:
     file("res_ref_summary.tsv") into refSummaryChannel
     file("res_tneg_summary.tsv") into tnegSummaryChannel

     shell:
     """
     # Reference
     best_model_ref=\$(tail -n +2 !{summary_ref} |head -1 |cut -s -f1|sed -e  "s/\\r//g"|sed "s/.pdb//g")
     dope_ref=\$(tail -n +2 !{summary_ref} |head -1 |cut -s -f3|sed -e  's/\\r//g')
     molpdf_ref=\$(tail -n +2 !{summary_ref} |head -1 |cut -s -f2|sed -e  's/\\r//g')
     normalized_dope_ref=\$(tail -n +2 !{summary_ref} |head -1 |cut -s -f4|sed -e  's/\\r//g')
     GA341_score_ref=\$(tail -n +2 !{summary_ref} |head -1 |cut -s -f5|sed -e  's/\\r//g')
     zscore_ref=\$(tail -n +2 !{prosa_ref} |head -1 |awk '{print \$2}'|sed -e  's/\\r//g'  )
     maxsub_ref=\$(tail -n +2 !{proq_ref} |head -1 |awk '{print \$2}'| sed -e  's/\\r//g' )
     lgscore_ref=\$(tail -n +2 !{proq_ref} |head -1 |awk '{print \$3}'| sed -e  's/\\r//g')
     echo -e "\$best_model_ref\t!{fam}\t\$molpdf_ref\t\$dope_ref\t\$normalized_dope_ref\t\$GA341_score_ref\t\$zscore_ref\t\$maxsub_ref\t\$lgscore_ref" > res_ref_summary.tsv
     # Negative
     best_model_tneg=\$(tail -n +2 !{summary_tneg} |head -1 |cut -s -f1|sed -e  "s/\\r//g"|sed "s/.pdb//g")
     dope_tneg=\$(tail -n +2 !{summary_tneg} |head -1 |cut -s -f3|sed -e  's/\\r//g')
     molpdf_tneg=\$(tail -n +2 !{summary_tneg} |head -1 |cut -s -f2|sed -e  's/\\r//g')
     normalized_dope_tneg=\$(tail -n +2 !{summary_tneg} |head -1 |cut -s -f4|sed -e  's/\\r//g')
     GA341_score_tneg=\$(tail -n +2 !{summary_tneg} |head -1 |cut -s -f5|sed -e  's/\\r//g')
     zscore_tneg=\$(tail -n +2 !{prosa_tneg} |head -1 |awk '{print \$2}'|sed -e  's/\\r//g'  )
     maxsub_tneg=\$(tail -n +2 !{proq_tneg} |head -1 |awk '{print \$2}'| sed -e  's/\\r//g' )
     lgscore_tneg=\$(tail -n +2 !{proq_tneg} |head -1 |awk '{print \$3}'| sed -e  's/\\r//g')
     echo -e "\$best_model_tneg\t!{fam}\t\$molpdf_tneg\t\$dope_tneg\t\$normalized_dope_tneg\t\$GA341_score_tneg\t\$zscore_tneg\t\$maxsub_tneg\t\$lgscore_tneg" > res_tneg_summary.tsv
     """
}

// run compute_structural_alignment.sh
process structural_alignment {
    tag "${fasta.baseName}:${fam}"
    publishDir "$myDir/modelling/${fam}_candidates/", mode: 'copy'
    cpus params.cpu

    input:
    set val(fam), file(fasta), best_pdb_ref, best_pdb_tneg from structuralAlignmentChannel

    output:
    file("*/struct_matrix_*.tsv") into structAliChannel
    file("*_candidates_ref_vs_ref/struct_matrix_mammoth.tsv") into refMammothChannel
    file("*_candidates_ref_vs_ref/struct_matrix_TMalign.tsv") into refTMalignChannel
    file("*_candidates_tneg_vs_tneg/struct_matrix_mammoth.tsv") into tnegMammothChannel
    file("*_candidates_tneg_vs_tneg/struct_matrix_TMalign.tsv") into tnegTMalignChannel

    shell:
    """
    mkdir !{fam}_candidates_ref_vs_ref/ !{fam}_candidates_tneg_vs_tneg/
    PDBRMSD.py -l !{fam} -q !{best_pdb_ref} -t ${params.database}/!{fam}/ref_pdb/ -s TMalign mammoth  -p !{params.tmalign_dir} !{params.mammoth_dir}  -r !{fam}_candidates_ref_vs_ref/ -b 1 -n !{params.cpu}
    PDBRMSD.py -l !{fam} -q !{best_pdb_tneg} -t ${params.database}/!{fam}/tneg_pdb/ -s TMalign mammoth  -p !{params.tmalign_dir} !{params.mammoth_dir}  -r !{fam}_candidates_tneg_vs_tneg/ -b 1 -n !{params.cpu}
    """
}


modelling_by_ref = refSummaryChannel.collectFile(name: 'candidates_by_ref.tsv')
modelling_by_ref.into { modelling_by_ref_tosave; modelling_by_ref_toprocess}
modelling_by_ref_tosave.subscribe { it.copyTo(modDir) }

modelling_by_tneg = tnegSummaryChannel.collectFile(name: 'candidates_by_tneg.tsv')
modelling_by_tneg.into { modelling_by_tneg_tosave; modelling_by_tneg_toprocess}
modelling_by_tneg_tosave.subscribe { it.copyTo(modDir) }

mammoth_by_ref = refMammothChannel.collectFile(name: 'mammoth_by_ref.tsv')
mammoth_by_ref.into { mammoth_by_ref_tosave; mammoth_by_ref_toprocess}
mammoth_by_ref_tosave.subscribe { it.copyTo(modDir) }

tmalign_by_ref = refTMalignChannel.collectFile(name: 'TMalign_by_ref.tsv')
tmalign_by_ref.into { tmalign_by_ref_tosave; tmalign_by_ref_toprocess}
tmalign_by_ref_tosave.subscribe { it.copyTo(modDir) }

mammoth_by_tneg = tnegMammothChannel.collectFile(name: 'mammoth_by_tneg.tsv')
mammoth_by_tneg.into { mammoth_by_tneg_tosave; mammoth_by_tneg_toprocess}
mammoth_by_tneg_tosave.subscribe { it.copyTo(modDir) }

tmalign_by_tneg = tnegTMalignChannel.collectFile(name: 'TMalign_by_tneg.tsv')
tmalign_by_tneg.into { tmalign_by_tneg_tosave; tmalign_by_tneg_toprocess}
tmalign_by_tneg_tosave.subscribe { it.copyTo(modDir) }

process build_matrix {
    publishDir "$myDir/", mode: 'copy'

    input:
    file(mamref) from mammoth_by_ref_toprocess
    file(mamtneg) from mammoth_by_tneg_toprocess
    file(tmref) from tmalign_by_ref_toprocess
    file(tmtneg) from tmalign_by_tneg_toprocess
    file(modref) from modelling_by_ref_toprocess
    file(modtneg) from modelling_by_tneg_toprocess

    output:
    file("pcm_result.tsv") into matrixChannel

    shell:
    """
    compute_matrix.py  -sr !{modref} -sn !{modtneg} -ar !{mamref} !{tmref} -an !{mamtneg} !{tmtneg} -o pcm_result.tsv
    """
}

process lineartest {
    publishDir "$myDir/", mode: 'copy'

    input:
    file(matrixChannel) from matrixChannel

    output:
    set file("reference_output.pdf"), file("ref_output.tsv"), file("prediction_output.tsv") into resultChannel

    shell:
    """
    Rscript !{baseDir}/bin/runlineartest2.R !{matrixChannel} !{params.universal_model} !{params.number_bootstrap} reference_output.pdf ref_output.tsv prediction_output.tsv
    """
}

println "Project : $workflow.projectDir"
println "Cmd line: $workflow.commandLine"
