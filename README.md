
# PCM : Pairwise Comparative Modelling
[Amine Ghozlane](https://research.pasteur.fr/fr/member/amine-ghozlane/) (amine.ghozlane@pasteur.fr) (@xealf8)

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Running PCM](#running-pcm)
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

For this example we will search resistance genes in the proteome of the following species:
<img src="example/phylogeny.png" align="center" />
This set was obtained from NCBI. The computation can be performed with the following command:
```
nextflow pcm.nf  --in example/example_proteome.faa --out result -w work/ -with-singularity pcm.img
```
Nextflow uses configuration file to deploy computation on cluster, an example of this file is available [here](nextflow_global.config).
This file enable with the profile singularity the usage of singularity on a slurm scheduler and need to be adjusted for each cluster configuration. Profiles are the activated with the command (-profile singularity  -c nextflow_global.config).

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
Ruppé E.¹, Ghozlane A.¹, Tap J.¹ et al. (2018) Prediction of the intestinal resistome by a novel 3D-based method. Nature Microbiology.
¹equally contributing authors.

## Contact

If you have any comments, questions or suggestions, or need help to use PCM, please contact the author [here](mailto:amine.ghozlane@pasteur.fr).
