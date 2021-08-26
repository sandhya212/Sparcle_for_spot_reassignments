# Sparcle for spot reassignments in transcriptomic images
This is the code base for Sparcle: Spatial reassignment of spots to cells via maximum likelihood estimation in transcriptomic images.

## Abstract

### Background
Imaging-based spatial transcriptomics has the power to reveal patterns of single-cell gene expression by detecting mRNA transcripts as individually resolved spots in multiplexed images. However, molecular quantification has been severely limited by the computational challenges of segmenting poorly outlined, overlapping cells, and of overcoming technical noise; the majority of transcripts are routinely discarded because they fall outside the segmentation boundaries. This lost information leads to less accurate gene count matrices and weakens downstream analyses, such as cell type or gene program identification. 
### Results
Here, we present Sparcle, a probabilistic model that reassigns transcripts to cells based on gene covariation patterns and incorporates spatial features such as distance to nucleus. We demonstrate its utility on both multiplexed error-robust fluorescence in situ hybridization (MERFISH) and single-molecule FISH (smFISH) data. 
### Conclusions
Sparcle improves transcript assignment, providing more realistic per-cell quantification of each gene, better delineation of cell boundaries, and improved cluster assignments. Critically, our approach does not require an accurate segmentation and is agnostic to technological platform.


## Graphical abstract
<img src=/figures/GUI_Mistic.png width="80%"></img>
## Datasets used 

1. MERFISH 
  - Merfish barcode and genes and paired scRNAseq data were obtained from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113576
  - FoVs were obtained from the authors (Moffitt et al) directly
  - Sparcle_ver_1.ipynb (or Sparcle_ver_1.py) is the parallelised Sparcle code to handle the Merfish FoVs
  - http://www.sciencemag.org/cgi/pmidlookup?view=long&pmid=30385464
  [Moffitt, J. R. Bambah-Mukku, D., Eichhorn, S. W., Vaughn, E., Shekhar, K., Perez, J. D. & Zhuang, X. Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. Science 2018.]


2. Allen smFISH VISp: 
 - https://zenodo.org/record/3478502#.YGSMi68zZ3g. 
 - Link courtesy from SSAM (J Park, W Choi, S Tiesmeyer, B Long, LE Borm, E Garren, TN Nguyen, S Codeluppi, M Schlesner, B Tasic, R Eils, N Ishaque. "Segmentation-free inference of cell types from in situ transcriptomics data." bioRxiv 800748. doi: https://doi.org/10.1101/800748)

3. STARmap: 
 - In this example we used the count matrix and tables for “visual_160/20171120_BF4_light" downloaded from: https://www.dropbox.com/sh/f7ebheru1lbz91s/AADiIArV5LmqzxdLvxo9qHXFa/visual_160/20171120_BF4_light?dl=0&subfolder_nav_tracking=1
 - Single-cell RNAseq data is from: Tasic B, et al. (2016) Adult mouse cortical cell taxonomy revealed by single cell transcriptomics.
(Nat Neurosci, 40, 335-346. DOI doi:10.1038/nn.4216) https://github.com/AllenInstitute/tasic2016data 

4. pciSeq/ISS: 
 - Data and visualisation code: https://colab.research.google.com/github/acycliq/pciSeq/blob/master/notebooks/pciSeq.ipynb
 - Slice used: CA1DapiBoundaries_4-3_right.tif from https://figshare.com/s/88a0fc8157aca0c6f0e8?file=13160426
 - Refer the pciSeq_ISS folder for how we process the data and run it through Sparcle

5. Vizgen FoV 75: 
 - We have used data from the recent ‘Vizgen Data release program’: https://info.vizgen.com/mouse-brain-data?utm_campaign=Data%20Release%20Program&utm_medium=email&_hsmi=125521312&_hsenc=p2ANqtz-_2jSArOjMVrw1OG1OC_o7XyhAuKxIq8oQ1d8NkKM8Cxn97U6P8rthj9kYCWVZ3JnGjMnuyXOvve7h42soJBmfnWEhrQxgJ3WOmhZ5W5E97l1JNki8&utm_content=125521312&utm_source=hs_email. 
- We have modified code provided in the corresponding Google colab to view the data for FoV 75

6. Code to extract the dangling mRNA for smFISH, STARmap, ISS and Vizgen will be made available upon request.

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

5. For another example of Sparcle, a standalone folder with data and code for pciSeq_ISS are provided.
