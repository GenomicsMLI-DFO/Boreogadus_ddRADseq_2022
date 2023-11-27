# Info --------------------------------------------------------------------

#  Check for sex-linked loci in Boreogadus saida
#  with the program RADsex
#
# Audrey Bourret
# 2023-11-17
#

# https://sexgenomicstoolkit.github.io/html/radsex/example.html

# Library -----------------------------------------------------------------

library(tidyverse)
library(Hmisc)
library(adegenet)

# Installation from GitHub requires devtools
#install.packages("devtools")
#devtools::install_github("SexGenomicsToolkit/sgtr")

devtools::install_github("biodray/sgtr")
library(sgtr)


# Data --------------------------------------------------------------------

pop.data <- read_csv("00_Data/00_FileInfos/Project_Infos_20230822.csv")
pop.data 


Arctogadus.ID <- c("S_22_00047", "S_22_00054", "S_22_00056", "S_22_00057", "S_22_00058", "S_22_00154", "S_22_00156",
          "S_22_00157", "S_22_00172")

load( file.path("./00_Data/06b_Filtering.ref", "B_10X_samples", "07_Final", "populations.38131snps.507indwArctogadus.adegenet.Rdata"))

gl.final


ID.Boreogadus <- pop.data %>% dplyr::filter(ID_GQ %nin% c(Arctogadus.ID ),
                                            ID_GQ %in% indNames(gl.final),
                                             Sexe_visuel %in% c("M", "F")) %>% 
  select(ID_GQ, Sexe_visuel)



ID.Boreogadus 

alig.files <- list.files("00_Data/03a_Demultiplex", full.names = T, pattern = "1.fq.gz")
alig.files




# Transfert files --------------------------------------------------------------

res.path <- "./02_Results/99_RADsex/Boreogadus/"

if(!dir.exists(res.path )){
  dir.create(res.path)
}

for(x in ID.Boreogadus$ID_GQ){
  print(x)
  ori.file <- alig.files %>% str_subset(pattern = paste0(x, ".1.fq.gz"))
  link.file <- paste0(res.path, x,".fq.gz")
  
  file.exists(link.file)
  
  cmd <- paste("-s", paste0("../../../",ori.file), link.file)
  print(cmd)
  system2("ln", cmd, stdout=T, stderr=T)
  
}

system2("ls", res.path)


# RADsex: Process ---------------------------------------------------------


cmd <- paste("process",
             "--input-dir",  res.path, 
             "--output-file", "./02_Results/99_RADsex/markers_table_Boreogadus.tsv",
             "--threads",  16)

cmd

A <- system2("radsex", cmd, stdout=T, stderr=T)
A

cat(file =  "./02_Results/99_RADsex/RADsex.process.log",
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")



cmd <- paste("depth", 
             "--markers-table",  "./02_Results/99_RADsex/markers_table_Boreogadus.tsv", 
             "--output-file", "./02_Results/99_RADsex/depth_Boreogadus.tsv",
             "--popmap", "./02_Results/99_RADsex/popmap_Boreogadus.tsv")
cmd

A <- system2("radsex", cmd, stdout=T, stderr=T)
A

cat(file =  "./02_Results/99_RADsex/RADsex.depth.log",
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")



sgtr::radsex_depth("./02_Results/99_RADsex/depth_Boreogadus.tsv", 
                     output_file = "./02_Results/99_RADsex/depth_Boreogadus.png")




# RADsex: Distrib ---------------------------------------------------------

# Creating a popmap

write.table(ID.Boreogadus, 
            file = "./02_Results/99_RADsex/popmap_Boreogadus.tsv",
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)


cmd <- paste("distrib", 
             "--markers-table",  "./02_Results/99_RADsex/markers_table_Boreogadus.tsv", 
             "--output-file", "./02_Results/99_RADsex/distribution_Boreogadus_10X.tsv",
             "--popmap", "./02_Results/99_RADsex/popmap_Boreogadus.tsv",
             "--min-depth", 10,
             "--groups M,F")
cmd

A <- system2("radsex", cmd, stdout=T, stderr=T)
A

cat(file =  "./02_Results/99_RADsex/RADsex.distrib_10X.log",
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")


sgtr::radsex_distrib("./02_Results/99_RADsex/distribution_Boreogadus_10X.tsv", 
                     output_file = "./02_Results/99_RADsex/distribution_Boreogadus_10X.png")

# RADsex: Signif ----------------------------------------------------------


cmd <- paste("signif",
             "--markers-table",  "./02_Results/99_RADsex/markers_table_Boreogadus.tsv", 
             "--output-file", "./02_Results/99_RADsex/significant_markers_table_Boreogadus_10X.tsv",
             "--popmap", "./02_Results/99_RADsex/popmap_Boreogadus.tsv",
             "--min-depth", 10,
             "--groups M,F"#, 
             #"--disable-correction"
)

A <- system2("radsex", cmd, stdout=T, stderr=T)
A


cat(file =  "./02_Results/99_RADsex/RADsex.signif.log",
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")


sgtr::radsex_marker_depths("./02_Results/99_RADsex/significant_markers_table_Boreogadus_10X.tsv", 
                           output_file = "./02_Results/99_RADsex/significant_markers_table_Boreogadus_10X.png", 
                           group_info_file = "./02_Results/99_RADsex/popmap_Boreogadus.tsv")


# Subset ------------------------------------------------------------------


test <- sgtr::load_marker_depths( "./02_Results/99_RADsex/significant_markers_table_Boreogadus_10X.tsv")


test %>% pivot_longer(cols = names(test)[3:ncol(test)], names_to = "ID_GQ", values_to = "N") %>% 
  dplyr::filter(N > 0) %>% 
  left_join(pop.data %>% dplyr::select(ID_GQ, Sexe_visuel, Region_echantillonnage, Longueur_mm, Numero_unique_groupe, Responsable_extraction)) %>% #head()
  #dplyr::filter(#Sexe_visuel == "F"
  #              Longueur_mm > 100
  #              ) %>% 
  ggplot(aes(x = ID_GQ, y = factor(id), fill =N)) +
  geom_bin2d() +
  scale_fill_viridis_c() +
  facet_grid(.~ Sexe_visuel +Responsable_extraction, scale = "free", space = "free") +
  theme(axis.text.x = element_blank())

              
ID.Boreogadus %>% left_join(pop.data) %>% 
  ggplot(aes(x = Longueur_mm)) + 
  geom_histogram() +
  facet_grid(Sexe_visuel ~ Region_echantillonnage)


cmd <- paste("subset",
             "--markers-table",  "./02_Results/99_RADsex/significant_markers_table_Boreogadus_5X.tsv", 
             "--output-file", "./02_Results/99_RADsex/subset_markers_table_Boreogadus_5X_females.tsv",
             "--popmap", "./02_Results/99_RADsex/popmap_Boreogadus.tsv",
             "--min-depth", 5,
             "--groups M,F", 
             "--min-group2", 5
             #"--disable-correction"
)

A <- system2("radsex", cmd, stdout=T, stderr=T)
A

sgtr::radsex_marker_depths("./02_Results/99_RADsex/subset_markers_table_Boreogadus_5X_females.tsv", 
                           output_file = "./02_Results/99_RADsex/significant_markers_table_Boreogadus_female.png", 
                           group_info_file = "./02_Results/99_RADsex/popmap_Boreogadus.tsv")




# RADsex: Map -------------------------------------------------------------

# Align to a chromosome

cmd <- paste("map", 
             "--markers-file", "./02_Results/99_RADsex/markers_table_Boreogadus.tsv",
             "--output-file", "./02_Results/99_RADsex/map_results_Boreogadus_5X.tsv",
             "--popmap", "./02_Results/99_RADsex/popmap_Boreogadus.tsv",
             "--genome-file", "./00_Data/99_REF_Genome/gadMor3.0/GCA_902167405.1/GCA_902167405.1_gadMor3.0_genomic.fna",
             "--min-depth", 5,
             "--groups M,F"#, 
)

cmd
A <- system2("radsex", cmd, stdout=T, stderr=T)
A

cat(file =  "./02_Results/99_RADsex/RADsex.map_5X.log",
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")

# Create a circo!!

radsex_map_circos("./02_Results/99_RADsex/map_results_Boreogadus_10X.tsv", 
                  detect_chromosomes = T,
                  chromosomes_file = "./02_Results/99_RADsex/gadus_chromosomes.tsv", 
                  output_file = "./02_Results/99_RADsex/map_results_Boreogadus_10X.png")


radsex_map_circos("./02_Results/99_RADsex/map_results_Boreogadus_5X.tsv", 
                  detect_chromosomes = T,
                  chromosomes_file = "./02_Results/99_RADsex/gadus_chromosomes.tsv", 
                  output_file = "./02_Results/99_RADsex/map_results_Boreogadus_5X.png")


radsex_map_region("./02_Results/99_RADsex/map_results_Boreogadus_10X.tsv", 
                  region = "chr13",
                  tracks = c("p", "bias"),
                  detect_chromosomes = T,
                  chromosomes_file = "./02_Results/99_RADsex/gadus_chromosomes.tsv", 
                  output_file = "./02_Results/99_RADsex/chr13_results_Boreogadus_10X.png"
           
)


radsex_map_region("./02_Results/99_RADsex/map_results_Boreogadus_5X.tsv", 
                  region = "chr13",
                  tracks = c("p", "bias"),
                  p_color = "#1E9B8A",
                  detect_chromosomes = T,
                  chromosomes_file = "./02_Results/99_RADsex/gadus_chromosomes.tsv", 
                  output_file = "./02_Results/99_RADsex/chr13_results_Boreogadus.pdf"
                  
)



plot("#FF0000")

radsex_map_manhattan( "./02_Results/99_RADsex/map_results_Boreogadus.tsv", 
                      #region = "chr13",
                      detect_chromosomes = T,
                      chromosomes_file = "./02_Results/99_RADsex/gadus_chromosomes.tsv", 
                      output_file = "./02_Results/99_RADsex/manhattan_results_Boreogadus.pdf")





if(!file.exists(file.path("./02_Results/99_RADsex/", ".gitignore")) ){
  cat("marker*.tsv", "*.pdf", "!.gitignore", "Boreogadus/*.fq.gz", sep = "\n",
      file = file.path("./02_Results/99_RADsex/", ".gitignore")) 
}





map.res <- sgtr::load_marker_depths( "./02_Results/99_RADsex/map_results_Boreogadus_5X.tsv")
head(map.res)


map.res %>% dplyr::filter(Signif == T) %>% arrange(Position) %>% head()
map.res %>% dplyr::filter(Signif == T) %>% arrange(Position) %>% tail() 


