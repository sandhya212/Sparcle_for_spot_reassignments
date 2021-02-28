# Sparcle for spot reassignments in transcriptomic images
This is the code base for Sparcle: Spatial reassignment of spots to cells via maximum likelihood estimation in transcriptomic images.

## Abstract

### Background
Imaging-based spatial transcriptomics has the power to reveal patterns of single-cell gene expression by detecting mRNA transcripts as individually resolved spots in multiplexed images. However, molecular quantification has been severely limited by the computational challenges of segmenting poorly outlined, overlapping cells, and of overcoming technical noise; the majority of transcripts are routinely discarded because they fall outside the segmentation boundaries. This lost information leads to less accurate gene count matrices and weakens downstream analyses, such as cell type or gene program identification. 
### Results
Here, we present Sparcle, a probabilistic model that reassigns transcripts to cells based on gene covariation patterns and incorporates spatial features such as distance to nucleus. We demonstrate its utility on both multiplexed error-robust fluorescence in situ hybridization (MERFISH) and single-molecule FISH (smFISH) data. 
### Conclusions
Sparcle improves transcript assignment, providing more realistic per-cell quantification of each gene, better delineation of cell boundaries, and improved cluster assignments. Critically, our approach does not require an accurate segmentation and is agnostic to technological platform.

## Datasets used
  * MERFISH data can be found here: http://www.sciencemag.org/cgi/pmidlookup?view=long&pmid=30385464
  [Moffitt, J. R. Bambah-Mukku, D., Eichhorn, S. W., Vaughn, E., Shekhar, K., Perez, J. D. & Zhuang, X. Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. Science 2018.]

  * Allen data: smFISH and matching scRNA-seq data of the primary visual cortex (VISp) region of an adult mouse brain can be found here: https://portal.brain-map.org/atlases-and-data/rnaseq#Mouse_Cortex_and_Hip. There are 3500 cells and 22 genes. 

## Installation

1. Download this code repository or Open Terminal and use `git clone`

 `$ git clone https://github.com/sandhya212/Sparcle_for_spot_reassignments`

2. The folder ‘Code_HPC’ contains the Python code implementing SPARCLE

* Sparcle_submit.sh: Shell script to submit code to the cluster
* start_file.py: Set current working directory and path variables for data here
* init_file.py: Initialise the variables for Sparcle here
* Sparcle.py: This file is Sparcle’s engine that iterates over FoVs recovering dangling mRNAs while parallel processing across FoVs per iteration

3. Submit code using: 

`
bsub -W 2:00 -R 'rusage[mem=25]' < Sparcle_submit.sh -o sparcle_out_file
`
4. For general reference, 'Code_HPC' equivalent is provided as a Python notebook and Python code in Sparcle_ver_1.ipynb and Sparcle_ver_1.py, respectively.
