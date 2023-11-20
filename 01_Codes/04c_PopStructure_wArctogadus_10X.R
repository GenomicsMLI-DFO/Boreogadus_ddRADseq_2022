# Info --------------------------------------------------------------------

# Population genetics analysis of Boreogadus saida
#  at global spatial scale
#  549 individuals (9 Arctogadus, 45 intermediate, 500 Boreogadus)
#
# Audrey Bourret
# 2023-11-17
#

# Library -----------------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)

library(adegenet)
library(hierfstat)

`%nin%` = Negate(`%in%`)

library(QuickPop)

# Data --------------------------------------------------------------------

pop.data <- read_csv("00_Data/00_FileInfos/Project_Infos_20230822.csv")
pop.data 

names(pop.data)

# Genetic Data ------------------------------------------------------------

load( file.path("./00_Data/06b_Filtering.ref", "B_10X_samples", "07_Final", "populations.38131snps.507indwArctogadus.adegenet.Rdata"))

gl.final
gi.final

# 
# na.gi.count <- function(gi){
#   res <- apply(tab(gi), MARGIN = 2, FUN = function(l){   n.na <- length(l[is.na(l) == T])
#   freq.na <- n.na / length(l)
#   return(freq.na)
#   })
#   res <- res[str_ends(names(res), "[.]0")] 
#   
#   names(res) <- names(res) %>% str_remove("[.]0")
#   
#   return(res)
#   
# }
# 
# # Function to create a list of loci, from a genind object
# 
# filter.MAF.NA <- function(gi, MAF.trs = 0.5, NA.trs = 0.5){
#   # Create vectors for each loci
#   MAF.res <- adegenet::minorAllele(gi)
#   NA.res  <- na.gi.count(gi)
#   
#   # Filter by threshold
#   MAF.loc <- dimnames(MAF.res[MAF.res >= MAF.trs])[[1]]
#   cat("There is", length( MAF.loc), "loci with MAF =", MAF.trs, "\n")
#   
#   NA.loc <- names(NA.res[NA.res <= NA.trs])
#   cat("There is", length(NA.loc), "loci with NA =", NA.trs, "\n")
#   
#   # LOCI with both conditions
#   LOCI.res <- c(MAF.loc, NA.loc)[duplicated(c(MAF.loc, NA.loc)) == T]
#   LOCI.res %>% length()
#   
#   cat("There is", length(LOCI.res), "loci with BOTH MAF =", MAF.trs, "and NA =" , NA.trs, "\n")
#   
#   return(LOCI.res)
# }
# 
# LOC.MAF05.NA10 <- filter.MAF.NA(gi.final, MAF.trs = 0.05, NA.trs = 0.10)
# LOC.MAF10.NA05 <- filter.MAF.NA(gi.final, MAF.trs = 0.10, NA.trs = 0.05)
# 
# 

# PCA ---------------------------------------------------------------------

pca.all  <- glPca(gl.final, center = TRUE, scale = FALSE,  
                   parallel = TRUE, n.core =16, nf = 1000)

#save(list = c("pca.all" ),
#     file = here("02_Results/01_PopStruct", "00_mtDNA", "PCA_2.Rdata"))

load( here::here("02_Results/01_PopStruct", "00_mtDNA", "PCA_2.Rdata"))

QuickPop::pca_varplot(pca.all)

Arctogadis.ID <- read_csv("02_Results/01_PopStruct/00_mtDNA/ArctogadusSeq1.csv")

pca.all %>% QuickPop::pca_scoretable(naxe = 10) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(Arctogadis.ID, by = c("ID" = "ID_GQ")) %>% 
  mutate(Taxon =  ifelse(is.na(Taxon), "Boreogadus", 
                  ifelse(ID %in% c("S_22_00047", "S_22_00054", "S_22_00056", "S_22_00057", "S_22_00058", "S_22_00154", "S_22_00156",
                                      "S_22_00157", "S_22_00172"), "Arctogadus", "Arcto-Boreo-gadus")),
         Region = ifelse(Region_echantillonnage %in% c("0A", "0B", "1A"), "0AB1A",
                  ifelse(Region_echantillonnage %in% c("2G", "2H", "2J"), "2GHJ", 
                         Region_echantillonnage ) 
                         )) %>% 
  dplyr::filter(Taxon !="Arctogadus") %>% 
  ggplot(aes(x = score.PC2, y = score.PC3, col = Region)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(~Region, ncol = 3) +
  stat_ellipse()+
  geom_point(alpha = 0.5, size = 2) + 
  #scale_color_viridis_d()+
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC2 (", QuickPop::pca_var(pca.all)$p.eig[2] %>% round(3) *100, "%)"),
    y = paste0("PC3 (", QuickPop::pca_var(pca.all)$p.eig[3] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) #+ theme(legend.position = "none")

