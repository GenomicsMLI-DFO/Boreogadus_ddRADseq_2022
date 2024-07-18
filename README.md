# R scripts and results of the Boreogadus genomic MS
## Genomic architecture and population structure of Boreogadus saida in Canadian waters: signals of adaptive divergence and hybridization with Arctogadus glacialis

__Main author:__  Trevor Bringloe and Audrey Bourret  
__Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        Laboratory of genomics   
__Location:__     Maurice Lamontagne Institute  
__Affiliated publication:__  Bringloe TT, Bourret A, Cote D, Roux MJ, Herbig J, Robert D, Geoffroy M, Parent GJ (accepted). enomic architecture and population structure of Boreogadus saida in Canadian waters: signals of adaptive divergence and hybridization with Arctogadus glacialis. Scientific Reports.
__Contact:__      e-mail: tbringloe@gmail.com

- [Abstract](#abstract)
- [Status](#status)
- [Contents](#contents)
- [Requirements](#requirements)
- [Acknowledgements](#acknowledgements)
- [References](#references)

## Abstract
The polar cod, Boreogadus saida, is an abundant and ubiquitous forage fish and a crucial link in Arctic marine trophic dynamics. Our objective was to unravel layers of genomic structure in B. saida from Canadian waters, specifically screening for potential hybridization with the Arctic cod, Arctogadus glacialis, large chromosomal inversions and sex-linked regions, prior to interpreting population structure. Our analysis of 53,384 SNPs in 522 individuals revealed hybridization and introgression between A. glacialis into B. saida. Subsequent population level analyses of B. saida using 12,305 SNPs in 511 individuals revealed featured three large (ca. 7.4 to 16.1 Mbp) chromosomal inversions, and a 2 Mbp region featuring sex-linked loci. We demonstrated population structuring across the Western and Eastern Arctic, and subarctic regions ranging from the Hudson Bay to the Atlantic maritime provinces. Genomic signal for the inferred population structure was highly aggregated into a handful of SNPs (13.8%), pointing to potentially important adaptive evolution across the Canadian range. Our study provides a high-resolution perspective on the genomic structure of B. saida, providing a foundation for work that could be expanded to the entire circumpolar range for the species.

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

Main script to obtain a first catalog of SNPs through STACKS. An older version of the script (02_Discovered_SNPs_2:30:10_old.R) using less longer reads was used in the first attempts to work on this dataset, but was later abandoned. 

#### 03c_Filter_SNPs_10X_2024.R

Main script to perform the first round of filtration, up to a SNPs panel with 1 SNP / locus, still with Arctogadus and hybrids individuals. Older version of the script (03a_Filter_SNPs_5X_2023.R and 03b_Filter_SNPs_10X_2023.R) were used in previous version of the pipeline.

#### 99_RADsex.R

Script to run the external program RADsex, aiming to identify sexualluy biased zone in the genome.  An older version of the script (99_RADsex_2:30:10_old.R) using less longer reads was used in the first attempts to work on this dataset, but was later abandoned. 

## Requirements
All sequence data can be accessed via NCBI’s BioProject PRJNA1062734. 

## Acknowledgement
Field work for this research was conducted aboard the CCGS Amundsen, the CCGS Teleost, the CCGS Leim, the Aqviq, the Katsheshuk II vessels and was supported by their crews as well as technical staff from the Marine Institute of Memorial University, Fisheries and Oceans Canada (DFO), Amundsen Science, the Nunatsiavut Government, the Université du Québec à Rimouski, and the Saguenay St-Lawrence Marine Park, Northern Shrimp Research Foundation. We thank Wojciech Walkusz (DFO) for organizing and managing the sample collection in the Hudson Strait region. The technical expertise for the laboratory work of Grégoire Cortial, Gabriel Bardaxoglou, Cloé Lepage, and Sara Khan (DFO) was also crucial valuable the success of this project.

## Reference





  
