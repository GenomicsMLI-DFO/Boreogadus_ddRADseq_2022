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
- Align on *GadMor3 (Gadus morhua)* reference genome. Did some tests with a B. saida assembly. The % alignment was lower (91% vs 98%) but in the end a similar number of SNPs was gentotyped (test on 26 individuals). 
  - Keep only those with > 96% alignment (574- remove 2 samples)
- Stacks v2.64 (stacks)
  - The catalog was made with the full dataset (574 samples)
  - **807,969** loci, effective per-sample coverage: mean=27.8x, stdev=14.7x, min=2.7x, max=94.8x
