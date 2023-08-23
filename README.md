# Boreogadus_ddRADseq_2022
ddRADseq on 500 Boreogadus


__Main author:__  Audrey Bourret  
__Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        Laboratory of genomics   
__Location:__     Maurice Lamontagne Institute  
__Affiliated publication:__  
__Contact:__      e-mail: audrey.bourret@dfo-mpo.gc.ca


- [Objective](#objective)
- [Summary](#summary)
- [Status](#status)
- [Contents](#contents)
  + [Subsections within contents](#subsections-within-contents)
- [Methods](#methods)
  + [Subsections within methods](#subsections-within-methods)
- [Main Results](#main-results)
- [Requirements](#requirements)
- [Caveats](#caveats)
- [Uncertainty](#uncertainty)
- [Acknowledgements](#acknowledgements)
- [References](#references)


## Objective
- Check population structure 

## Summary
Description of the project, provide some background and context. What are the inputs and outputs?


## Status
In-development


## Contents
### Folder structure

```
.
├── 00_Data                 # Main datasets are here (large files are stored locally)
├── 01_Code                 # R Scripts for the different analysis
├── 02_Results              # Figures and main output results
└── README.md
```

### Main scripts 

#### 02_Discovered_SNPs.R

Main script to obtain a first catalog of SNPs through STACKS


### Subsections within contents
Use subsections to describe the purpose of each script if warranted


## Methods
### SNP panels

**Main catalog**

- 6 librairies of 96 samples = 576, sequenced on 1 lane of NovaSeq 150 PE
- Remove Illumina adaptor with Trimmomatic  
- Demultiplex with process_radtag module, 2 restriction enzymes (*pstI* & *mspI*)), truncation at 135 pb, filter-illumina, clean, rescue and check quality options
- Align on *GadMor3 (Gadus morhua)* reference genome. Did some tests with a B. saida assembly. The % alignment was lower (91% vs 98%) but in the end a similar number of SNPs was genotyped (test on 26 individuals). 
  - Keep only those with > 96% alignment (574- remove 2 samples)
- Stacks v2.64 (stacks)
  - The catalog was made with the full dataset (574 samples)
  - **807,969** loci, effective per-sample coverage: mean=27.8x, stdev=14.7x, min=2.7x, max=94.8x

**High quality SNPs panel**

  - Keep only samples with mean min coverage = 5x (562 samples)
  - Populations (stacks)
      - Parameters: R = 0.75 overall, min MAF = 0.01
      - Kept 54,247 loci and 211,357 snps
  - Remove missing values (individuals and SNPS)
      - VCFtools - 0.1.17
      - Few individuals with more than 30% missing values, remove duplicated individuals, choosing the replicate with the less missing data 
      - X SNPs with more than 10% missing values
      - After filtration : 170,307 snps from 553 individuals
  - HW desequilibrium by sampling location - not performed
      - Instead, remove X loci with He > 0.6 
  - Each SNPs depth - Remove a few SNPs with too low or high median coverage (< x and >x, 2 times SD, approx. 1-99% percentile)
      - Remove  SNPs (x snps from x individuals)
  - Check relatedness
      - With VCFtool, relatedness2
      - Remove x pairs of samples with relatedness > 0.25 (~0.45, so clearly duplicates)
      - final : 3398 individuals
      - 9 individuals with relatedness << 1 , also outliers within a PCA. Maybe Arctogadus glacialis.
  - Keep only 1 snps by RADloc (Highest MAF / lowest NA), MAF 0.01 (some were fixed)
      - **28,309 snps from 540 individuals**


