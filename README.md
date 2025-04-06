# AD-KSIF (Assembly-Driven KeystoneSpecies Identification Framework)

## 1. Overview of AD-KSIF

Development of keystone species identificationoffers fundamental force for transformative synthetic microbial community applicationto meet the controllability demand of modern engineering systems, yet the lackof underlying mechanisms support and limited information for functionalcategorization build barriers in module embedding. Raising the key idea ofhighlighting keystone species by characteristic assembly trajectories, thenallocating the species into annotated functional categories, here we propose anovel assembly-driven keystone species identification framework (AD-KSIF) toresolve the challenge. Validated by real community datasets from 6 watertreatment engineering systems, AD-KSIF outperforming the mainstream baselinemodel and exhibits robustness against the community stochasticity. The outputsfrom AD-KSIF enable us to pick cross-system conserved keystone species againstthe disruption of insufficient domestication, and to link them into diversefunctional modules of minimal genome-scale metabolic reconstruction. Thepresented framework demonstrates the sensitivity of environment-communityinteraction marks and the power of introducing transparent principles intonumeric ecological modelling, paving the way for the systematically data-drivenmanagement of next-generation engineering synthetic microbial community design.

## 2. Preparation for AD-KSIF pipeline

### 2.1 R packages

AD-KSIF pipeline is developed in R language (v 4.4.2), which relies on the packages as follow:

- ape

- dplyr

- doParallel

- foreach

- muscle

- picante

- stringr

Please ensure that you have installed these packages correctly before using the pipeline, and we appreciately integrated the script from co-authur Qilin Wang for phylogenetic binning module in function source 'PRE_trt.R'. The original code is available on [GitHub - Kylin-Wang/trefile2binningOTU: Bin the OTUs in the tre file according to phylogenetic distance](https://github.com/Kylin-Wang/trefile2binningOTU/tree/main).

### 2.2 Directory and data preparation

AD-KSIF pipeline allows using raw datasets of OTU abundance table (in .csv format) and sequences (in .fasta format) for complete keystone species identification. In order to make it easier for users to adjust, We set amount of intermediate outputs carrying information from each module of the pipeline. Therefore, we suggest a well-organized directory heirachy for storation. Here is an example:

/AD-KSIF data

    |__/Functions

        |__PARA_cal.R

        |__PRE_trt.R

        |__TRE_crt.R

    |__/Input_data

        |__OTU-abundance.csv

        |__OTU_seq.fasta

    |__Bins

    |__Indexes

    |__Results

Then we are ready for an exploration by AD-KSIF.

## 3. Introduction of AD-KSIF pipeline modules

### 3.1 Functions

In order to make the code simple and easy to understand, we modularize the pipeline and encapsulate some of the core function in R source files. Here's a brief introduction of them.

1. TRE_crt.R.
   This file contains the functions for phylogenetic tree construction.
   *Tree_Construction(df): the function accept data frame read by ‘readDNAStringSet' in 'muscle' package, and do multi-sequence alignment in algorithm of muscle and phylogenetic tree construction in NJ method. The phylogenetic tree is for return.

2. PRE_trt.R.
   This file contains the functions for binning.
   *fasta_treatment(df): the function accept data frame of FASTA file and divide the ID and the sequence to prepare suitable format of input data. A two-colomn dataframe of 'ID' and 'Sequence' is for return.
   
   *binning(d.max,n.min,input): the function accept data frame return form 'fasta_treatment' and do binning according to the given parameters of phylogenetic distance threshold (d.max) and the minimum size of a bin (n.min). A two-colomn dataframe of 'ID' and 'Bin' is for return.
   
   *allocation(otubin,otuseq,otuabund): the function combines the OTUs, OTU sequences, OTU bins and OTU reads in each sample for subsequent analysis.
   
   The detail information of these functions is available on [GitHub - Kylin-Wang/trefile2binningOTU: Bin the OTUs in the tre file according to phylogenetic distance](https://github.com/Kylin-Wang/trefile2binningOTU/tree/main). 

3. PARA_cal.R.
   This file contains the functions for index calculation.
   
   *beta_NTI(otu,phylo,n): the function accept OTU abundance table and phylogenetic tree for beta_NTI calculation in each bin at the iteration number of n.
   
   *RC_BC(otu,phylo,n): the function accept OTU abundance table for tRC_Bray-Curtis calculation in each bin at the iteration number of n.
   
   *X_Y(otu,bNTI,RC): the function accept OTU abundance table, beta_NTI table and RC_Bray-Curtis table in each bin for evolution and dispersal index calculation.
   
   *beta_NTI(otu,phylo,n): the function accept evolution and dispersal index table in each bin for stability and category index calculation.

### 3.2 Main pipeline

The main pipeline is divided into 6 modules with outputs from each:

- load the function files

- data pre-treatment

- phylogenetic tree construction

- binning

- index calculation

- result output

## 4. Implementation of AD-KSIF pipeline

It's easy to run the pipeline with the script 'AD-KSIF_pipeline.R' just by R studio, bash, or any tools you like. Just remember to fill the blank path for each dictionary and the target samples you are interested in the whole table. You can find the blanks just in the script. Have a good time with it!

## 5. Citation

The work is still under review and not published now. Any form of citation need a contact with Kong Boning with the email adress of kbn23@mails.tsinghua.edu.cn .


