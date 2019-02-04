
# PCM : Pairwise Comparative Modelling
[Amine Ghozlane](https://research.pasteur.fr/fr/member/amine-ghozlane/) (amine.ghozlane@pasteur.fr) (@xealf8)

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Running PCM](#running-pcm)
- [Results](#results)
- [Dependencies](#dependencies)
- [Citation](#citation)
- [Contact](#contact)


## Introduction

PCM is a generic method using homology modelling to increase the specificity of functional prediction of proteins, especially
when they are distantly related from proteins for which a function is known. The principle of PCM is to
build structural models and assess their relevance using a specific training approach. PCM uses the
list of sequences of reference proteins from a given family, the structures related to this family (they
will be used as structural templates in the PDB format) and a series of negative references.
The pcm process follow this dag:
<img src="flowchart.png" align="center" />


## Installation

PCM is deployed using [singularity](https://singularity.lbl.gov/) and nextflow.
To install singularity on a Linux platform, follow these commands:
```
VERSION=2.5.2
wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
tar xvf singularity-$VERSION.tar.gz
cd singularity-$VERSION
./configure --prefix=/usr/local
make
sudo make install
```
Download PCM singularity image (warning: the file is heavy ~ 19Go):
```
wget ftp://shiny01.hosting.pasteur.fr/pub/pcm.img
```
Install nextflow:
```
curl -s https://get.nextflow.io | bash
```

## Running PCM

### Start

For this example we will search resistance genes in the proteome of the following species:
<img src="example/phylogeny.png" align="center" />
This set was obtained from NCBI. The computation can be performed with the following command:
```
git clone https://github.com/aghozlane/pcm.git
cd pcm
# See help with
nextflow pcm.nf --help
# Run example calculation
nextflow pcm.nf  --in example/example_proteome.faa --out result -w work/ -with-singularity pcm.img
```

The fasta file must be in the following format:
```
>ID1
MNTFGQIHNNMPYLFLLAFIMNFYDQFNNSISGQEMCYEVESI
FNNHQVDIIGAPAAAFKPLELQKGLGTKGAIVNYPILQVTGNI
>ID2
MNTFGQIHNNMPYLFLLAFIMNFYDQFNNSISGQEMCYEVESI
FNNHQVDIIGAPAAAFKPLELQKGLGTKGAIVNYPILQVTGNI
```
The ID is >name (without any space).

### Singularity

The usage of singularity with Nextflow requires that the input data/result/work to be accessible in the virtual image. If you do not set your calculation in your home folder, you might need to explain to nextflow where it can mount this space, like here:
```
#nextflow config file 
singularity {                                                                    
            enabled = true                                                       
            autoMounts = false                                                   
            runOptions = '--bind  /mnt:/home'                                     
}
# the directory /mnt that contains my data is mounted in the singularity 
# on its home directory

# Do not run
# You can then provide this configuration like this:
nextflow pcm.nf  --in example/example_proteome.faa --out result -w work/ -with-singularity pcm.img -c nextflow.config
```
More information about nextflow and singularity are available [here](https://www.nextflow.io/docs/latest/singularity.html)

### Cluster configuration

Nextflow uses a configuration file to deploy computation on cluster, an example of this file is available [here](nextflow_global.config). 
This file enables the usage of singularity on a slurm scheduler and need to be adjusted for each cluster configuration. Profiles are the activated with the command (-profile singularity  -c nextflow_global.config).
```
# The profile singularity correspond to an execution on a cluster environment
profiles{
    singularity {
        process.container = 'file:///pasteur/homes/aghozlan/pcm/img/pcm.img'
        singularity.enabled = true
        singularity.autoMounts = false
        singularity.runOptions = "-B /pasteur:/pasteur -B /local:/local -H ~/:/mnt"
        process.executor = 'slurm'
        process.clusterOptions='--qos=hubbioit'
        process.queue = 'hubbioit'
    }
}
# Do not run
# You can then provide this configuration like this:
nextflow pcm.nf  --in example/example_proteome.faa --out result -w work/ -with-singularity pcm.img -c nextflow_global.config -profile singularity
```
In the profile singularity, I indicate where is the container (with process.container). I mount data directory with singularity.runOptions (the path are here indicatives). I select the type of scheduler (with executor) and our cluster option: the qos and the partition for hubbioit.

## Results

PCM results are available as following:

File | Directory | Description
---|---|---
**pcm_result.tsv** | main directory | This table gives all the score obtained during modelisation with reference and negative reference modelling. This file is provided to our classifier to determine if the candidate is more a reference or a negative reference.
**prediction_output.tsv** | main directory | This table is the result of the classification of each candidates. It provides the pvalue associated to the score to be a potential reference and the estimated quality of the prediction.
**ref_output.tsv** | main directory |  This table is the result of the classification for the universal model sequence (same as prediction_output.tsv).
**reference_output.pdf** | main directory | This pdf shows the roc curve and weight associated to each score during the learning phase of machine learning. It uses an "universal model" (/usr/local/bin/database/universal_model.csv) that contains the score that we obtained for know ARD and know "not ARD".
**_candidates.fasta** | candidates | This fasta file  gives the selected hit from blast/ssearch/hmmer and length criterion for one ARD family
**_candidates_hit_length.tsv** | candidates | It indicates which hit meet the length criterion 
**_candidates_hit_properties** | candidates | It regroups all the information of blast/ssearch/hmmer and clustalo alignments for each candidates (Identity,Similarity, coverage)
**candidates_by_ref.tsv** | modelling | This table gives the results obtained for candidates modelled with a reference template in the following order: molpdf,	dope_score, normalized_dope,	GA341_score,	zscore,	maxsub_ref,	lgscore
**candidates_by_tneg.tsv** | modelling | This table gives the results obtained for candidates modelled with a negative reference template (same as ref models)
**mammoth_by_ref.tsv** | modelling | This table gives of the structural alignment results with the reference templates using Mammoth (candidate ard_type, PDB_reference,	zscore_mammoth,	TM-score_mammoth)
**mammoth_by_tneg.tsv** | modelling | This table gives the structural alignment results with the negative reference templates using Mammoth (same as mammoth by ref)
**TMalign_by_ref.tsv** | modelling | This table gives of the structural alignment results with the reference templates using TMalign (candidate ard_type,	PDB_reference,	RMSD_TMalign,	TM-score_TMalign)
**TMalign_by_tneg.tsv** | modelling | This table gives the structural alignment results with the negative reference templates using TMalign (same as TMalign by ref)
**best_model_ref** | modelling/*_candidates/ | This directory provide the best model obtained for one given candidate using a reference template
**best_model_tneg** | modelling/*_candidates/ | This directory provide the best model obtained for one given candidate using a negative reference template

The other files correspond to intermediate results.

## Dependencies

*  __Blastp__
  Search of distant homologous for the candidate sequence.
   _Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. "Basic local alignment search tool." J. Mol. Biol., 1990, 215:403-410._
* __Clustal 0mega__
Perform multiple alignment between candidate and reference protein sequence.
	Sievers, F., & Higgins, D. G. (2014). Clustal Omega, accurate alignment of very large numbers of sequences. In _Multiple sequence alignment methods_ (pp. 105-116). Humana Press, Totowa, NJ.
* __HMMER__
 Search of distant homologous for the candidate sequence.
Johnson, L. S., Eddy, S. R., & Portugaly, E. (2010). Hidden Markov model speed heuristic and iterative HMM search procedure. _BMC bioinformatics_, _11_(1), 431.
* __Modeller__
 Performs homology modelling of candidate sequence and reference protein structures.
 Fiser, A., & Do, R. K. G. (2000). Modeling of loops in protein structures. _Protein science_, _9_(9), 1753-1773.
* __ProQ__
Performs quality checking protein models.
Wallner, B., & Elofsson, A. (2003). Can correct protein models be identified?. _Protein science_, _12_(5), 1073-1086.
* __Prosa__
Performs quality checking protein models.
Wiederstein, M., & Sippl, M. J. (2007). ProSA-web: interactive web service for the recognition of errors in three-dimensional structures of proteins. _Nucleic acids research_, _35_(suppl_2), W407-W410.
* __Psipred__
Predicts the secondary structure of candidate sequence. This secondary secondary structure is used by modeller as a topology constraint.
Jones, D. T. (1999). Protein secondary structure prediction based on position-specific scoring matrices1. _Journal of molecular biology_, _292_(2), 195-202.
* __ssearch__
Search of distant homologous for the candidate sequence.
Pearson, W. R., & Lipman, D. J. (1988). Improved tools for biological sequence comparison. _Proceedings of the National Academy of Sciences_, _85_(8), 2444-2448.
* __R packages__
	* __LiblinearR__
Performs logist regression to automatically classify protein models.
Helleputte, T. LiblineaR: Linear Predictive Models Based on the LIBLINEAR C/C++ Library, 2015. _R package version_, 1-94.
	* __pROC__
Draw a classification ROC curves
Robin, X., Turck, N., Hainard, A., Tiberti, N., Lisacek, F., Sanchez, J. C., & Müller, M. (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. _BMC bioinformatics_, _12_(1), 77.

## Citation

Please cite the following publication if you use PCM:
Ruppé E.¹, Ghozlane A.¹, Tap J.¹ et al. (2018) Prediction of the intestinal resistome by a novel 3D-based method. Nature Microbiology (Article available [here](https://www.nature.com/articles/s41564-018-0292-6)).  


## Contact

If you have any comments, questions or suggestions, or need help to use PCM, please contact the author [here](mailto:amine.ghozlane@pasteur.fr).
