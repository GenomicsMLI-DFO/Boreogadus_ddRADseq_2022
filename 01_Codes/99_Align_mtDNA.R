# Info --------------------------------------------------------------------

# Small scrip to align a few samples to B saida mtDNA
#
# Audrey Bourret
# 2023-08-18
#

# Library -----------------------------------------------------------------

library(parallel)
library(tidyverse)
library(readxl)


# Internal functions
for(i in 1:length( list.files("./01_Codes/Functions") )){
  source(file.path("./01_Codes/Functions",  list.files("./01_Codes/Functions")[i]))  
}

`%nin%` = Negate(`%in%`)

# Add pythom env to this specific project
Sys.setenv(PATH = paste(c("/home/genyoda/Documents/PythonVenv/GenoBaseEnv/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))

system2("multiqc", "--help")

# Check folders -----------------------------------------------------------
# To run if you need the folder skeleton 

#auto.folder()

# Data --------------------------------------------------------------------

pop.info <- read.table(file.path(get.value("info.path"), "popmap.txt"))
names(pop.info) <- c("Sample", "POP")
pop.info


pop.data <- read_csv(file.path(get.value("info.path"),"Project_Infos_20230811.csv"))

pop.data 


# Define initial working directory (just in case something goes wrong)
current.wd <- getwd()

numCores <- if(get_os() %in% c("os","linux")){
  detectCores() # Utilise le max de coeurs si sur linux
} else 1

cat("There's", numCores, "cores available in this computer", sep = " ")

# BWA - index the reference genome ----------------------------------------

# Try first with an old Boreogadus assembly
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_900302515.1/

# Check that 
A <- system2("bwa", "", stdout=T, stderr=T)

list.files(file.path(get.value("ref.genome.path"), "mtDNA_Bsaida"))

cmd <- paste("index",
             "-p",  file.path(get.value("ref.genome.path"),"mtDNA_Bsaida"), 
             "-a", "bwtsw",
             file.path(get.value("ref.genome.path"),"mtDNA_Bsaida", "NC_010121.1.fasta")
)


A <- system2("bwa", cmd, stdout=T, stderr=T)
A

# save a log file 

cat(file = file.path(get.value("ref.genome.path"), "Bsaida.mtDNA.August2023.log" ),
    cmd, "\n\n",
    A, # what to put in my file
    append= T, sep = "\n")



# Perform the alignment ---------------------------------------------------



# Complete version
demulti.files <- list.files(get.value("demulti.path"), pattern = ".1.fq.gz", full.names = T) %>% 
  str_subset(".rem.", negate = T)# %>%  
 # str_subset("S_22_00047|S_22_00054|S_22_00056|S_22_00057|S_22_00058|S_22_00154|S_22_00156|S_22_00172|S_22_07689|S_22_09427|S_22_00230|S_22_00089")

demulti.files %>% length()

#demulti.files[1:10]

mclapply(demulti.files,
         FUN = function(x){
           # How I will rename all this : 
           file.R1 <- x
           file.R2 <- x %>% str_replace(".1.fq.gz", ".2.fq.gz")
           file.bam <- x %>% str_replace(get.value("demulti.path"),  "./00_Data/03c_Align_mtDNA") %>% 
             str_replace(".1.fq.gz", ".bam")
           file.sort.bam <- file.bam %>% str_replace(".bam", ".sorted.bam")
           stat.tsv <- file.bam %>% str_replace(".bam", ".stat.tsv")
           
           # DO THE ALIGMENT  
           if(!file.exists(stat.tsv)){ # do it only if it doesn't exist
             
             cmd1 <- paste("mem",
                           "-t", 16,
                           "-M",
                           
                           file.path(get.value("ref.genome.path"),"mtDNA_Bsaida"), # the index ref genome
                           file.R1,
                           file.R2,
                           "2> /dev/null",
                           "| samtools", "view", "-F 4", 
                           "--threads", 15, # Number of additional threads
                           "-o", file.bam#,
             )# the file
             
             A <- system2("bwa", cmd1, stdout=T, stderr=T)
             
             # Sort the file 
             cmd2 <- paste("sort",
                           "--threads", 15,
                           #"-n",
                           "-O", "BAM",
                           # the new file (sorted)
                           "-o", file.sort.bam,
                           # the bam file to sort
                           file.bam
                           
             )
             
             A <- system2("samtools", cmd2, stdout=T, stderr=T)
             
             
             # cmd2.5 <- paste("view",
             #                 file.sort.bam,
             #               "| head -n 100 >",
             #                 file.sort.bam %>% str_replace(".sorted.bam", ".subset.sam")
             #               # the bam file to sort
             #               
             #               
             # )
             # 
             # A <- system2("samtools", cmd2.5, stdout=T, stderr=T)
             # 
             
             # Compute stats
             
             cmd3 <- paste("flagstat",
                           "--threads", 15,
                           "-O", "tsv",
                           file.sort.bam,
                           ">",
                           stat.tsv
             )
             
             A <- system2("samtools", cmd3, stdout=T, stderr=T)
             
           }
           
         },
         mc.cores = 2
) 



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")

library("Rsamtools")


which <- GRanges(c(
  "NC_010121.1:1000-2000"
))
## equivalent:
## GRanges(
##     seqnames = c("seq1", "seq2", "seq2"),
##     ranges = IRanges(
##         start = c(1000, 100, 1000),
##         end = c(2000, 1000, 2000)
##     )
## )
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)


bam.files <- list.files("00_Data/03c_Align_mtDNA", 
                        pattern = ".sorted.bam",
                        full.names = T)


seq1 <-  DNAStringSet()
seq2 <-  DNAStringSet()
seq3 <-  DNAStringSet()

for(x in 1:length(bam.files)){

  print(bam.files[x])
  
  bam.int <- scanBam(bam.files[x])    
  
  seq1.int <- bam.int[[1]]$seq[which(bam.int[[1]]$pos == 4894)]
  seq2.int <- bam.int[[1]]$seq[which(bam.int[[1]]$pos == 5138)]
  seq3.int <- bam.int[[1]]$seq[which(bam.int[[1]]$pos == 9134)]
  
  if(length(seq1.int) > 0) {
    
  seq1.cons.int <- consensusString(seq1.int) 
  seq1 <- c(seq1, seq1.cons.int)  
  names(seq1)[length(seq1)] <- bam.files[x] %>% str_remove("00_Data/03c_Align_mtDNA/") %>% str_remove(".sorted.bam")  
  }
  
  if(length(seq2.int) > 0) {
    seq2.cons.int <- consensusString(seq2.int)
    seq2 <- c(seq2, seq2.cons.int)  
    names(seq2)[length(seq2)] <- bam.files[x] %>% str_remove("00_Data/03c_Align_mtDNA/") %>% str_remove(".sorted.bam") 
  }
  
  if(length(seq3.int) > 0) {
    seq3.cons.int <- consensusString(seq3.int)
    seq3 <- c(seq3, seq3.cons.int)  
    names(seq3)[length(seq3)] <- bam.files[x] %>% str_remove("00_Data/03c_Align_mtDNA/") %>% str_remove(".sorted.bam") 
  }
  
  
  }


seq1
seq2
seq3



library(Biostrings)

writeXStringSet(seq1, filepath = "Seq1_mtDNA.fasta")
writeXStringSet(seq2, filepath = "Seq2_mtDNA.fasta")
writeXStringSet(seq3, filepath = "Seq3_mtDNA.fasta")



# Blast -------------------------------------------------------------------


# Load taxonomy data

ncbi.tax <- readr::read_tsv("/media/genyoda/Fast_Storage/Projet/Banque_REF_NCBI.git/01_Data/01_blastdb/rankedlineage.dmp", 
                            col_names = c("id", "name", "species", "genus", "family", "order", "class","phylum", "kingdom", "superkingdom"), 
                            col_types=("i-c-c-c-c-c-c-c-c-c-"))


NCBI.path <- "/media/genyoda/Fast_Storage/Projet/Banque_REF_NCBI.git/01_Data/01_blastdb/"

cmd1 <- paste("-db", "nt", 
              "-query",  file.path("Seq1_mtDNA.fasta"),
              "-outfmt", "\"7 qseqid sacc staxid ssciname sskingdom pident length mismatch gapopen qstart qend sstart send evalue bitscore\"",
              "-out", file.path("Blast.Seq1.out"), 
              "-perc_identity", 95,
              "-num_threads", 8,
              #"-max_target_seqs", 10, 
              sep = " ")# forward adapter

#A<-system2("blastn", cmd1, stdout=T, stderr=T,
#           env = paste0("BLASTDB=", NCBI.path))
A


RES.Seq1.ncbi <- read.table(file.path("Blast.Seq1.out"),
                           sep="\t")


names(RES.Seq1.ncbi) <- c("QueryAccVer", "SubjectAccVer", "TaxoId","SciName", "SKindom", "Identity", "AlignmentLength", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")

RES.Seq1.ncbi %>% head()

RES.Seq1.ncbi <- RES.Seq1.ncbi %>% left_join(ncbi.tax, by = c("TaxoId" = "id"))


cmd1 <- paste("-db", "nt", 
              "-query",  file.path("Seq2_mtDNA.fasta"),
              "-outfmt", "\"7 qseqid sacc staxid ssciname sskingdom pident length mismatch gapopen qstart qend sstart send evalue bitscore\"",
              "-out", file.path("Blast.Seq2.out"), 
              "-perc_identity", 95,
              "-num_threads", 8,
              #"-max_target_seqs", 10, 
              sep = " ")# forward adapter

#A<-system2("blastn", cmd1, stdout=T, stderr=T,
#           env = paste0("BLASTDB=", NCBI.path))
A


RES.Seq2.ncbi <- read.table(file.path("Blast.Seq2.out"),
                            sep="\t")


names(RES.Seq2.ncbi) <- c("QueryAccVer", "SubjectAccVer", "TaxoId","SciName", "SKindom", "Identity", "AlignmentLength", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")

RES.Seq2.ncbi <- RES.Seq2.ncbi %>% left_join(ncbi.tax, by = c("TaxoId" = "id"))

cmd1 <- paste("-db", "nt", 
              "-query",  file.path("Seq3_mtDNA.fasta"),
              "-outfmt", "\"7 qseqid sacc staxid ssciname sskingdom pident length mismatch gapopen qstart qend sstart send evalue bitscore\"",
              "-out", file.path("Blast.Seq3.out"), 
              "-perc_identity", 95,
              "-num_threads", 8,
              #"-max_target_seqs", 10, 
              sep = " ")# forward adapter

#A<-system2("blastn", cmd1, stdout=T, stderr=T,
#           env = paste0("BLASTDB=", NCBI.path))
A


RES.Seq3.ncbi <- read.table(file.path("Blast.Seq3.out"),
                            sep="\t")


names(RES.Seq3.ncbi) <- c("QueryAccVer", "SubjectAccVer", "TaxoId","SciName", "SKindom", "Identity", "AlignmentLength", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")

RES.Seq3.ncbi <- RES.Seq3.ncbi %>% left_join(ncbi.tax, by = c("TaxoId" = "id"))




BLAST_LCA <- function(RES, threshold = 0.97){
  DF <- data.frame()
  
  RES.OK <- RES %>% filter(AlignmentLength >= .95 * width,
                           Identity >= threshold) %>% 
    mutate(Taxon = NA,
           Levels = NA)
  
  # if from ncbi
  if(str_count(names(RES.OK), "SciName") %>% sum() == 1){
    RES.OK <- RES.OK %>% filter(str_detect(SciName, "environmental sample|uncultured|predicted", negate = T))
    RES.OK$species <- paste(sapply(str_split(RES.OK$SciName, " "),`[`,1),
                            sapply(str_split(RES.OK$SciName, " "),`[`,2))
    RES.OK$genus <- sapply(str_split(RES.OK$specie, " "),`[`,1)
    RES.OK <- RES.OK %>% mutate(species = ifelse(str_detect(species, " sp[.]| cf[.]| aff."), NA, species))
  }
  
  ASV <- RES.OK %>% pull(QueryAccVer) %>% unique()
  for(x in seq_along(ASV)){
    RES.INT <- RES.OK %>% filter(QueryAccVer == ASV[x])
    
    # loops around ranks
    for(y in c("species", "genus", "family", "order", "class", "phylum", "kingdom")) {
      
      N.LCA <- RES.INT %>% filter(!is.na(y)) %>%  pull(y) %>% unique() %>% str_subset("NA", negate = T) %>% length()
      
      if(N.LCA==1){
        RES.INT$Taxon <-  RES.INT %>% filter(!is.na(y)) %>% pull(y) %>% unique()  %>% str_subset("NA", negate = T)
        RES.INT$Levels <-  y
        break;
      }else{
        RES.INT[,y] <- NA
      }
    } # END of the loop around columns
    #RES.INT <- RES.INT %>% select(QueryAccVer, kingdom, phylum, class, order, family, genus, species, Taxon, Levels) %>% 
    #                     distinct(.keep_all = T)
    
    if(nrow(RES.INT[which(!is.na(RES.INT[,y])),])>0){
      DF <- bind_rows(DF, RES.INT[which(!is.na(RES.INT[,y])),])      
    }
    
    
  }
  return(DF)  
  
}   

BLAST_TOPHIT <- function(RES, threshold = 0.95){
  DF <- data.frame()
  
  RES.OK <- RES %>% filter(#AlignmentLength >= .95 * width,
    Identity >= threshold) %>% 
    mutate(Taxon = NA,
           Levels = NA)
  
  # if from ncbi
  if(str_count(names(RES.OK), "SciName") %>% sum() == 1){
    RES.OK <- RES.OK %>% filter(str_detect(SciName, "environmental sample|uncultured|predicted", negate = T))
    RES.OK$species <- paste(sapply(str_split(RES.OK$SciName, " "),`[`,1),
                            sapply(str_split(RES.OK$SciName, " "),`[`,2))
    RES.OK$genus <- sapply(str_split(RES.OK$species, " "),`[`,1)
    RES.OK <- RES.OK %>% mutate(species = ifelse(str_detect(species, " sp[.]| cf[.]| aff."), NA, species))
  }
  
  ASV <- RES.OK %>% pull(QueryAccVer) %>% unique()
  for(x in seq_along(ASV)){
    RES.INT <- RES.OK %>% filter(QueryAccVer == ASV[x])
    
    evalue.min <- min(RES.INT$evalue)
    RES.INT <- RES.INT %>% filter(evalue == evalue.min)
    
    #identity.max <- max(RES.INT$Identity)
    #RES.INT <- RES.INT %>% filter(Identity == identity.max)
    
    # loops around ranks
    for(y in c("species", "genus", "family", "order", "class", "phylum", "kingdom")) {
      
      N.LCA <- RES.INT %>% filter(!is.na(y)) %>%  pull(y) %>% unique() %>% str_subset("NA", negate = T) %>% length()
      
      if(N.LCA==1){
        RES.INT$Taxon <-  RES.INT %>% filter(!is.na(y)) %>% pull(y) %>% unique()  %>% str_subset("NA", negate = T)
        RES.INT$Levels <-  y
        break;
      }else{
        RES.INT[,y] <- NA
      }
    } # END of the loop around columns
    #RES.INT <- RES.INT %>% select(QueryAccVer, kingdom, phylum, class, order, family, genus, species, Taxon, Levels) %>% 
    #                     distinct(.keep_all = T)
    
    if(nrow(RES.INT[which(!is.na(RES.INT[,y])),])>0){
      DF <- bind_rows(DF, RES.INT[which(!is.na(RES.INT[,y])),])      
    }
    
  }
  return(DF)  
  
}   


sum.BLAST <- function(DF){
  RES <- DF %>% select(QueryAccVer, Taxon, Levels, species, genus, family, order, class, phylum, kingdom) %>% unique() 
  
  return(RES)
}


Seq1.BlastTOP.95 <- BLAST_TOPHIT(RES.Seq1.ncbi, threshold = 99)  %>% sum.BLAST() %>% 
                          dplyr::select(ID_GQ =  QueryAccVer, Taxon ) %>%  mutate(Loc = "Seq1")
Seq2.BlastTOP.95 <- BLAST_TOPHIT(RES.Seq2.ncbi, threshold = 99)  %>% sum.BLAST() %>% dplyr::select(ID_GQ =  QueryAccVer, Taxon ) %>% mutate(Loc = "Seq2")
Seq3.BlastTOP.95 <- BLAST_TOPHIT(RES.Seq3.ncbi, threshold = 99)  %>% sum.BLAST() %>%dplyr::select(ID_GQ =  QueryAccVer, Taxon ) %>%  mutate(Loc = "Seq3")


head(Seq3.BlastTOP.95)

bind_rows(Seq1.BlastTOP.95,
          Seq2.BlastTOP.95,
          Seq3.BlastTOP.95) %>% #pivot_wider(names_from = Loc, values_from = Taxon) %>% 
  left_join(pop.data) %>% group_by(Loc, Region_echantillonnage, Taxon) %>% 
  summarise(N = n()
            ) %>% 
  mutate(Ntot = sum(N),
         Nprop = N/Ntot) %>%# head()
  ggplot(aes(x = Region_echantillonnage, y = Taxon, fill = Nprop)) +
  geom_bin2d() +
  scale_fill_viridis_c() +
  facet_grid(. ~ Loc) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pop.data %>% 
  left_join(Seq1.BlastTOP.95) %>% 
  ggplot(aes(x = Longueur_mm, fill = Taxon)) +
  geom_histogram() +
  #scale_fill_viridis_c() +
  facet_grid(Region_echantillonnage ~ .) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



bind_rows(Seq1.BlastTOP.95,
          Seq2.BlastTOP.95,
          Seq3.BlastTOP.95) %>% #pivot_wider(names_from = Loc, values_from = Taxon) %>% 
  left_join(pop.data) %>% #group_by(Loc, Region_echantillonnage, Taxon) %>% 
  #dplyr::filter(Region_echantillonnage == "Beaufort_Sea_ecoregion") %>% 
  #dplyr::mutate(CAT = ifelse(ID_GQ %in% c("S_22_00047","S_22_00054","S_22_00056","S_22_00057","S_22_00058","S_22_00154","S_22_00156","S_22_00172"), "Weird", "Other")) %>% 
  dplyr::mutate(Region_echantillonnage= ifelse(ID_GQ %in% c("S_22_00047","S_22_00054","S_22_00056","S_22_00057","S_22_00058","S_22_00154","S_22_00156","S_22_00172"), "Weird_nuclear", Region_echantillonnage)) %>% 
  
  ggplot(aes(x = Loc, y = ID_GQ, fill = Taxon)) +
  geom_bin2d() +
  #scale_fill_viridis_c() +
  facet_wrap(~Region_echantillonnage, scale = "free_y") +
  #facet_grid(CAT ~ ., scale = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


Arctogadis.ID <-  Seq1.BlastTOP.95 %>% dplyr::filter(Taxon == "Arctogadus glacialis")

write_csv(Arctogadis.ID, "./02_Results/01_PopStruct/00_mtDNA/ArctogadusSeq1.csv")

bind_rows(Seq1.BlastTOP.95,
          Seq2.BlastTOP.95,
          Seq3.BlastTOP.95) %>% #pivot_wider(names_from = Loc, values_from = Taxon) %>% 
  left_join(pop.data) %>% #group_by(Loc, Region_echantillonnage, Taxon) %>% 
  dplyr::filter(Region_echantillonnage %in% c("2H", "2J", "3K")) %>% 
  dplyr::mutate(CAT = ifelse(ID_GQ %in% c("S_22_00047","S_22_00054","S_22_00056","S_22_00057","S_22_00058","S_22_00154","S_22_00156","S_22_00172"), "Weird", "Other")) %>% 
  ggplot(aes(x = Loc, y = ID_GQ, fill = Taxon)) +
  geom_bin2d() +
  #scale_fill_viridis_c() +
  facet_grid(CAT ~ ., scale = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



Seq1.BlastTOP.95 %>% View()


RES.Seq1.ncbi %>% View()
