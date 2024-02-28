# Info --------------------------------------------------------------------

# PCADAPT for Boreogadus saida
#  at global spatial scale
#  549 individuals (9 Arctogadus, 45 intermediate, 500 Boreogadus)
#
# Audrey Bourret
# 2024-02-27
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

# Genetic Data Conversion ------------------------------------------------------------

vcf.path <- "/media/genyoda/Extra_Storage/Projets/Data_Trevor/02_Results/08_Boreogadus/New_analysis_23ii24/populations.11233snps.511indwArctogadus.H06.DP.single.final.recode.vcf"

vcf.data <- vcfR::read.vcfR(vcf.path)
gl.data  <- vcfR::vcfR2genlight(vcf.data) 
gi.data  <- vcfR::vcfR2genind(vcf.data) 

SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, ","), `[`,1) %>% str_remove("scaffold"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )

SCAFFOLD.info



#save(list = c("vcf.data", "gl.data", "gi.data"),
#     file ="TrevorVCF_2024.Rdata"
#     )

cmd1 <- paste("--vcf", vcf.path, 
              #"--recode",
              "--plink-tped",
              "--out",  "02_Results/02_PCADAPT_2024/populations.11233snps.511indwArctogadus.H06.DP.single.final.recode")


cmd1

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)
A1

cmd2a <- paste("--tfam", "./02_Results/02_PCADAPT_2024/populations.11233snps.511indwArctogadus.H06.DP.single.final.recode.tfam", 
               "--tped", "./02_Results/02_PCADAPT_2024/populations.11233snps.511indwArctogadus.H06.DP.single.final.recode.tped", 
               "--make-bed", 
               "--out", "./02_Results/02_PCADAPT_2024/populations.11233snps.511indwArctogadus.H06.DP.single.final" 
               
)

 A2a <- system2("/home/genyoda/Documents/Programs/plink_linux_x86_64_20210606/plink", cmd2a, stdout=T, stderr=T)
A2a


# PCADAPT -----------------------------------------------------------------

library(pcadapt)
#BiocManager::install("qvalue")
library("qvalue")

# Convertion to plink .bed format

# Read .bed in PCAadapt
pcadapt.genotype  <- read.pcadapt(file.path("./02_Results/02_PCADAPT_2024/populations.11233snps.511indwArctogadus.H06.DP.single.final.bed" ),
                                         type = "bed")


pcadapt.snp <- read.delim(file.path(here::here(), "./02_Results/02_PCADAPT_2024/populations.11233snps.511indwArctogadus.H06.DP.single.final.bim" ),
                                 header = F) %>% pull(V2)

# Run pcadapt

K.init <- 10

pcadapt.k10   <- pcadapt(pcadapt.genotype, K =K.init)
# Check screeplot

plot(pcadapt.k10, option = "screeplot") 

# Check structure

plot(pcadapt.k10, option = "scores") 

pcadapt.k2 <- pcadapt(pcadapt.genotype , K = 2)
pcadapt.k3 <- pcadapt(pcadapt.genotype , K = 3)
pcadapt.k4 <- pcadapt(pcadapt.genotype , K = 4)


plot(pcadapt.k2, option = "manhattan")

hist(pcadapt.k2$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

pcadapt.k2.qqplot <- plot(pcadapt.k2, option = "qqplot") + theme_bw()

pcadapt.k2.chi2 <- plot(pcadapt.k2, option = "stat.distribution")   + theme_bw()


ggsave(filename = here::here("02_Results/02_PCADAPT_2024/pcadapt_k2_chi2.pdf"), 
       plot =pcadapt.k2.chi2, 
       width =5, height =5 , units = "in",
       dpi = 300)

ggsave(filename = here::here("02_Results/02_PCADAPT_2024/pcadapt_k2_qqplot.pdf"), 
       plot =pcadapt.k2.qqplot, 
       width =5, height =5 , units = "in",
       dpi = 300)

# Statistics
#x$pvalues 
alpha <- 0.05 

qval.k2 <- qvalue::qvalue(pcadapt.k2$pvalues)$qvalues
outliers.k2 <- which(qval.k2 < alpha)
length(outliers.k2)

length(outliers.k2) / nLoc(gl.data)

padj.k2 <- p.adjust(pcadapt.k2$pvalues,method="BH")
outliers.padj.k2 <- which(padj.k2 < alpha)
length(outliers.padj.k2)

bonf.k2 <- p.adjust(pcadapt.k2$pvalues,method="bonferroni")
outliers.bonf.k2 <- which(bonf.k2 < alpha)
length(outliers.bonf.k2)

length(outliers.bonf.k2)/ nLoc(gl.data)

library(eulerr)
library(RColorBrewer)

venn.pcadapt <- plot(euler(list(qvalue =   pcadapt.snp[outliers.k2]  , BH = pcadapt.snp[outliers.padj.k2], Bonferroni =  pcadapt.snp[outliers.bonf.k2] ) , shape = "circle"),
                     quantities = T,
                     fill = brewer.pal(n = 3, name = "Dark2"),
                     alpha = 0.7, edges = F)

venn.pcadapt

res.pcadapt <- data.frame(ID = pcadapt.snp,
                          pvalues = pcadapt.k2$pvalues,
                          qvalues = qval.k2) %>% 
               left_join(SCAFFOLD.info %>% dplyr::select(ID, CHROM, POS)) %>% 
  mutate(Outliers.qvalue = ifelse(ID %in% pcadapt.snp[outliers.k2], T, F),
         Outliers.bonferroni = ifelse(ID %in% pcadapt.snp[outliers.bonf.k2], T, F))

res.pcadapt %>% View()

res.pcadapt %>% mutate(test = -log10(pvalues)) %>% arrange(pvalues) %>% View()

#write_csv(res.pcadapt, "02_Results/02_PCADAPT_2024/PCADAPT_k2_results.csv")

qval.k3 <- qvalue::qvalue(pcadapt.k3$pvalues)$qvalues
outliers.k3 <- which(qval.k3 < alpha)
length(outliers.k3)

res.pcadapt %>% dplyr::filter(str_detect(CHROM, "LR")) %>% 
                                group_by(CHROM) %>% 
                                summarise(Nsnps = n(),
                                          Noutlier = length(ID[Outliers.bonferroni == T]),
                                          Nneutral = length(ID[Outliers.bonferroni == F])) %>% 
  
                              mutate(Prop = Noutlier / Nsnps) %>% arrange(desc(Prop)) %>% 
  View()


qval.k4 <- qvalue::qvalue(pcadapt.k4$pvalues)$qvalues
outliers.k4 <- which(qval.k4 < alpha)
length(outliers.k4)
#outliers.high.sfa.k2 <- c(which(is.na(qval.sfa.k2)), which(qval.sfa.k2 < 3.5e-6)) 


pcadapt.outliers <- pcadapt.snp[outliers.k2]

#save(list = c("pcadapt.outliers"),
#     file = file.path("02_Results/02_PopulationGenetics/OutliersLoc.Rdata"))

nLoc(gl.data)


venn.pcadapt <- plot(euler(list(pcadapt.k2=   pcadapt.snp[outliers.k2]  , pcadapt.k3 =  pcadapt.snp[outliers.k3], pcadapt.k4 =  pcadapt.snp[outliers.k4] ) , shape = "circle"),
                     quantities = T,
                     fill = brewer.pal(n = 3, name = "Dark2"),
                     alpha = 0.7, edges = F)

venn.pcadapt


pcadapt.outliers <- pcadapt.snp[outliers.k2]


# Create new VCF ----------------------------------------------------------

pcadapt.outliers %>% length() / nLoc(gl.data)


df.outlier <-  SCAFFOLD.info %>% 
  dplyr::filter(ID %in% pcadapt.outliers)

df.neutral <-  SCAFFOLD.info %>% 
  dplyr::filter(ID %nin% pcadapt.outliers)

nrow(df.outlier)
nrow(df.neutral)


write.csv(df.outlier %>% select(ID), file.path("02_Results/02_PCADAPT_2024/Outlier_pcadapt_k2.csv"), 
          row.names = F, quote = F)

write.csv(df.neutral %>% select(ID), file.path("02_Results/02_PCADAPT_2024/Neutral_pcadapt_k2.csv"), 
          row.names = F, quote = F)

# CREATE VCF WITH UNIQUE

vcf.path 

cmd <- paste("--vcf", vcf.path, 
             "--recode",
             "--snps", file.path("02_Results/02_PCADAPT_2024/Outlier_pcadapt_k2.csv"),
             
             "--out", file.path("02_Results/02_PCADAPT_2024/populations.outlier")
             )


cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)


cmd <- paste("--vcf", vcf.path, 
             "--recode",
             "--snps", file.path("02_Results/02_PCADAPT_2024/Neutral_pcadapt_k2.csv"),
             
             "--out", file.path("02_Results/02_PCADAPT_2024/populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.noinver.neutral")
)


cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)


# PCA ---------------------------------------------------------------------

library(Hmisc)

pca.outliers  <- glPca(gl.data[,locNames(gl.data) %in%  pcadapt.snp[c(outliers.bonf.k2) %>% unique()]], center = TRUE, scale = FALSE,  
                  parallel = TRUE, n.core =16, nf = 1000)

pca.neutral  <- glPca(gl.data[,locNames(gl.data) %nin% pcadapt.outliers], center = TRUE, scale = FALSE,  
                       parallel = TRUE, n.core =16, nf = 1000)


pca.eco  <- glPca(gl.data[,locNames(gl.data) %in% bayes.ecoregion$LOC[bayes.ecoregion$qval < 0.05]], center = TRUE, scale = FALSE,  
                      parallel = TRUE, n.core =16, nf = 1000)


pca.outliers %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate( Region = ifelse(Region_echantillonnage %in% c("0A", "0B", "1A"), "0AB1A",
                                  ifelse(Region_echantillonnage %in% c("2G", "2H", "2J"), "2GHJ", 
                                         Region_echantillonnage ) 
  )) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Region)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(~Region_echantillonnage, ncol = 3) +
  stat_ellipse(aes(col = Region))+
  geom_point(alpha = 0.5, size = 2) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.outliers)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.outliers)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8)


pca.eco %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate( Region = ifelse(Region_echantillonnage %in% c("0A", "0B", "1A"), "0AB1A",
                          ifelse(Region_echantillonnage %in% c("2G", "2H", "2J"), "2GHJ", 
                                 Region_echantillonnage ) 
  )) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Region)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(~Region_echantillonnage, ncol = 3) +
  stat_ellipse(aes(col = Region))+
  geom_point(alpha = 0.5, size = 2) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.eco)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.eco)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8)



pca.neutral %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate( Region = ifelse(Region_echantillonnage %in% c("0A", "0B", "1A"), "0AB1A",
                          ifelse(Region_echantillonnage %in% c("2G", "2H", "2J"), "2GHJ", 
                                 Region_echantillonnage ) 
  )) %>% 
  #left_join(Env.data) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Region)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(~Region_echantillonnage, ncol = 3) +
  #stat_ellipse(aes(col = Region))+
  geom_point(alpha = 0.5, size = 2) +  
#  scale_color_viridis_c()+
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.neutral)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.neutral)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8)



# DIFFERENT PEAKS BY CLUSTER? ---------------------------------------------

# Step 1 - categorized them

clust.test <- snapclust.choose.k(x = gi.data,
                                     max = 6)
which.min(clust.test)

plot(clust.test, type = "b", cex = 2, xlab = "k", ylab = "AIC")
points(which.min(clust.test), min(clust.test), col = "blue", pch = 20, cex = 2)
abline(v = which.min(clust.test), lty = 2, col = "red")

clust.k3 <- snapclust(x = gi.data, k = 3)
clust.k4 <- snapclust(x = gi.data, k = 4)

k3.df <- data.frame(ID_GQ =  names(clust.k3$group),  
                    Clust =  clust.k3$group)

k4.df <- data.frame(ID_GQ =  names(clust.k4$group),  
                    Clust =  clust.k4$group)


# Step 1 - PCA on all

pca.all  <- glPca(gl.data, center = TRUE, scale = FALSE,  
                       parallel = TRUE, n.core =16, nf = 1000)


pca.all %>% QuickPop::pca_scoretable(naxe = 10) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(k4.df, by = c("ID" = "ID_GQ")) %>% 
  left_join(ID_Ecoregion, by = c("ID" = "ID_GQ")) %>% 
 ggplot(aes(x = score.PC1, y = score.PC2, col = Clust)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~Region_echantillonnage, ncol = 3) +
  stat_ellipse()+
  geom_point(alpha = 0.5, size = 2) + 
   labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.all)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.all)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) #+ theme(legend.position = "none")

# Step 2 - BayeScan 

Env.data <- pop.data %>%   
  group_by(Numero_unique_groupe) %>% 
  summarise(Lat = mean(Latitude_echantillonnage_DD),
            Long = mean(Longitude_echantillonnage_DD))

eco.shp <-  terra::vect("00_Data/00_FileInfos/SIG/FederalMarineBioregions_SHP/FederalMarineBioregions.shp")
env.shp <-  terra::vect(Env.data, geom = c("Long", "Lat"), keepgeom = T, crs = "+proj=longlat" )

terra::crs(eco.shp, proj = T)

ev2 <- terra::project(env.shp , terra::crs(eco.shp, proj = T))

terra::vect()

plot(eco.shp, "NAME_E")
plot(ev2, add = T, col = "darkgray")

res <- terra::intersect(ev2, eco.shp)

ecoregion.df <- as.data.frame(res) %>%  dplyr::select(Numero_unique_groupe,  NAME_E)


pop.data %>% left_join(ecoregion.df) %>% group_by(Region_echantillonnage, NAME_E) %>% summarise(N = n())

ID_Ecoregion <- pop.data %>% left_join(ecoregion.df) %>% dplyr::mutate(Ecoregion = ifelse(Region_echantillonnage == "1A", "Eastern Arctic",
                                                                   ifelse(Region_echantillonnage == "Saguenay", "Gulf of Saint Lawrence",
                                                                   ifelse(NAME_E == "Arctic Basin", "Western Arctic", NAME_E)))) %>% 
                                                         dplyr::select(ID_GQ, Ecoregion)

pop.data %>% left_join(ecoregion.df) %>% left_join(ID_Ecoregion) %>%   group_by(Region_echantillonnage, NAME_E, Ecoregion) %>% summarise(N = n())

#write_csv(ID_Ecoregion, "ID_ecoregion.csv")

# Step 3 - Bayscan
# Outlier can be loaded without reanalyses

library(dartR)
?gl2bayescan

# Data conversion
gl.data@other$loc.metrics.flags$monomorphs <- FALSE
gl.data@ploidy <- rep(as.integer(2), nInd(gl.data))


# Regional Bayescan

pop(gl.data) <- data.frame(ID_GQ = indNames(gl.data)) %>% 
  left_join(ID_Ecoregion) %>% pull(Ecoregion)

pop(gl.data) %>% table()

gl2bayescan(gl.data, outfile = "02_Results/03_BayeScan_2024/bayescan_gl.ecoregion.txt", outpath = ".")

# Clust 1 vs 2 vs 3

#ID.c1c2 <- k3.df %>% dplyr::filter(Clust %in% c(1,2)) %>% pull(ID_GQ)
#ID.c1c3 <- k3.df %>% dplyr::filter(Clust %in% c(1,3)) %>% pull(ID_GQ)
#ID.c2c3 <- k3.df %>% dplyr::filter(Clust %in% c(2,3)) %>% pull(ID_GQ)

#gl.c1c2 <- gl.data[indNames(gl.data) %in% ID.c1c2,]

#pop(gl.c1c2) <- data.frame(ID_GQ = indNames(gl.c1c2)) %>% 
#  left_join(k3.df) %>% pull(Clust)

#pop(gl.c1c2) %>% table()

#gl2bayescan(gl.c1c2, outfile = "02_Results/03_BayeScan/bayescan_gl.c1c2.txt", outpath = ".")

#gl.c1c3 <- gl.data[indNames(gl.data) %in% ID.c1c3,]

#pop(gl.c1c3) <- data.frame(ID_GQ = indNames(gl.c1c3)) %>% 
#  left_join(k3.df) %>% pull(Clust)

#pop(gl.c1c3) %>% table()

#gl2bayescan(gl.c1c3, outfile = "02_Results/03_BayeScan/bayescan_gl.c1c3.txt", outpath = ".")

#gl.c2c3 <- gl.data[indNames(gl.data) %in% ID.c2c3,]

#pop(gl.c2c3) <- data.frame(ID_GQ = indNames(gl.c2c3)) %>% 
#  left_join(k3.df) %>% pull(Clust)

#pop(gl.c2c3) %>% table()

#gl2bayescan(gl.c2c3, outfile = "02_Results/03_BayeScan/bayescan_gl.c2c3.txt", outpath = ".")

# Run bayescan

cmd <- paste("02_Results/03_BayeScan_2024/bayescan_gl.ecoregion.txt",
             "-snp",
             "-od", "./02_Results/03_BayeScan_2024/Results/", #Output directory file
             "-threads", 20
)

cmd

A <- system2("bayescan", cmd, stdout = T, stderr = T)
A
 
# cmd <- paste("02_Results/03_BayeScan/bayescan_gl.c1c2.txt",
#              "-snp",
#              "-od", "./02_Results/03_BayeScan/Results/", #Output directory file
#              "-threads", 20
# )
# 
# cmd
# 
# A <- system2("bayescan", cmd, stdout = T, stderr = T)
# A
# 
# cmd <- paste("02_Results/03_BayeScan/bayescan_gl.c1c3.txt",
#              "-snp",
#              "-od", "./02_Results/03_BayeScan/Results/", #Output directory file
#              "-threads", 20
# )
# 
# cmd
# 
# A <- system2("bayescan", cmd, stdout = T, stderr = T)
# A
# 
# cmd <- paste("02_Results/03_BayeScan/bayescan_gl.c2c3.txt",
#              "-snp",
#              "-od", "./02_Results/03_BayeScan/Results/", #Output directory file
#              "-threads", 20
# )
# 
# cmd
# 
# A <- system2("bayescan", cmd, stdout = T, stderr = T)
# A
# 


source("/home/genyoda/Documents/Programs/BayeScan/R_functions/plot_R.r")

plot_bayescan("02_Results/03_BayeScan_2024/Results/bayescan_gl.ecoregion_fst.txt", FDR = 0.05, add_text= F)
#plot_bayescan("02_Results/03_BayeScan/Results/bayescan_gl.c1c2_fst.txt", FDR = 0.05, add_text= F)
#plot_bayescan("02_Results/03_BayeScan/Results/bayescan_gl.c1c3_fst.txt", FDR = 0.05, add_text= F)
#plot_bayescan("02_Results/03_BayeScan/Results/bayescan_gl.c2c3_fst.txt", FDR = 0.05, add_text= F)


bayes.ecoregion <- read.table("02_Results/03_BayeScan_2024/Results/bayescan_gl.ecoregion_fst.txt")
#bayes.c1c2      <- read.table("02_Results/03_BayeScan/Results/bayescan_gl.c1c2_fst.txt")
#bayes.c1c3      <- read.table("02_Results/03_BayeScan/Results/bayescan_gl.c1c3_fst.txt")
#bayes.c2c3      <- read.table("02_Results/03_BayeScan/Results/bayescan_gl.c2c3_fst.txt")

add.bayes.names <- function(fst, gl){
  loc <- data.frame(LOC = adegenet::locNames(gl)) %>% 
    mutate(LOC.num = row.names(.))
  
  fst %>% mutate(LOC.num = row.names(.)) %>% left_join(loc)
  
}

hist(bayes.ecoregion$fst)
#hist(bayes.c1c2$fst)
#hist(bayes.c1c3$fst)
#hist(bayes.c2c3$fst)

# The estimated alpha coefficient indicating the strength and direction of selection. 
# A positive value of alpha suggests diversifying selection, whereas negative values 
# suggest balancing or purifying selection.
#

bayes.ecoregion <- add.bayes.names(bayes.ecoregion, gl.data)
#bayes.c1c2 <- add.bayes.names(bayes.c1c2, gl.data)
#bayes.c1c3 <- add.bayes.names(bayes.c1c3, gl.data)
#bayes.c2c3 <- add.bayes.names(bayes.c2c3, gl.data)

bayes.ecoregion %>% filter(qval < 0.05) %>%   pull(alpha) %>% hist() 

# How many SNPs identified by ecoregion

bayes.ecoregion %>% filter(qval <= 0.05) %>% nrow / bayes.ecoregion %>% nrow()
bayes.ecoregion %>% filter(qval <= 0.05) %>% nrow 

library("ggVennDiagram")


library(eulerr)


# Default plot
venn.outlier <- plot(euler(list(#ALL = locNames(gl.final),
  ecoregion =  bayes.ecoregion$LOC[bayes.ecoregion$qval < 0.05], 
  #c1c2 = bayes.c1c2$LOC[bayes.c1c2$qval < 0.05],
  #c1c3 =  bayes.c1c3$LOC[bayes.c1c3$qval < 0.05],
  #c2c3 = bayes.c2c3$LOC[bayes.c2c3$qval < 0.05],
 qvalue =   pcadapt.snp[outliers.k2]  #, 
 # Bonferroni =  pcadapt.snp[outliers.bonf.k2] 
) , shape = "circle"),
quantities = T,
#fill = brewer.pal(n = 3, name = "Dark2"),
alpha = 0.7, edges = T)

venn.outlier

venn.bayescan

446 / 11233
236/314

# Default plot
library("ggVennDiagram")
ggVennDiagram(list(#ALL = locNames(gl.final),
  BayeScan_Region =  bayes.ecoregion$LOC[bayes.ecoregion$qval < 0.05], 
  #c1c2 = bayes.c1c2$LOC[bayes.c1c2$qval < 0.05],
  #c1c3 =  bayes.c1c3$LOC[bayes.c1c3$qval < 0.05],
  #c2c3 = bayes.c2c3$LOC[bayes.c2c3$qval < 0.05],
  pcadapt_qvalue =   pcadapt.snp[outliers.k2]  , 
  pcadapt_bonferroni =  pcadapt.snp[outliers.bonf.k2] 
), label_alpha = 0, label = c("count")) +
   scale_fill_distiller(palette = "RdBu")

# 
# ggVennDiagram(list(#ALL = locNames(gl.final),
#   #BayeScan_Region =  bayes.ecoregion$LOC[bayes.ecoregion$qval < 0.05], 
#   c1c2 = bayes.c1c2$LOC[bayes.c1c2$qval < 0.05],
#   c1c3 =  bayes.c1c3$LOC[bayes.c1c3$qval < 0.05],
#   c2c3 = bayes.c2c3$LOC[bayes.c2c3$qval < 0.05]#,
#   #pcadapt_qvalue =   pcadapt.snp[outliers.k2]  , 
#   #pcadapt_bonferroni =  pcadapt.snp[outliers.bonf.k2] 
# ), label_alpha = 0, label = c("count")) +
#   scale_fill_distiller(palette = "RdBu")
# 
# 
# test <-   bind_rows( bayes.c1c2 %>% mutate(Comp = "c1c2"),
#   bayes.c1c3 %>% mutate(Comp = "c1c3"),
#   bayes.c2c3 %>% mutate(Comp = "c2c3"),
#   bayes.ecoregion %>% mutate(Comp = "Region")
# )%>% left_join(SCAFFOLD.info, by = c("LOC" = "ID")) %>% 
#   dplyr::mutate(CHROM = ifelse(!str_detect(CHROM, "LR"), "Unplaced", CHROM),
#                 POS = as.numeric(POS))  
# 
# test <- test %>% left_join(Chromo.df)
# test  %>% 
#   mutate(sign = ifelse(qval < 0.05/8000, "highly yes",
#                 ifelse(qval < 0.05, "yes", "no"))) %>% 
#   ggplot(aes(x = POS, y = fst, col = sign )) +
#   geom_point() +
#   facet_grid(Comp~chr, scale = "free_x", space = "free_x")+
#   theme_minimal() + 
#   theme(axis.text.x = element_blank(),
#         strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
#         strip.text.y = element_text(angle = 90),
#         panel.grid = element_blank(),
#         panel.spacing = unit(0, "cm"),
#         panel.border = element_rect(fill = NA, colour = "black"),
#         plot.background = element_rect(fill = "white", colour  = "white"),
#         legend.title = element_blank(),
#         axis.title.x = element_blank(),
#         legend.position = "bottom",
#         plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )
# 


library(qqman)

head(test)

Chromo.df <- data.frame(CHROM = unique(test$CHROM))
Chromo.df$chr <- 1:24


test %>% str()
# Make the Manhattan plot on the gwasResults dataset
manhattan(test , chr="chr", bp="POS", snp="LOC", p="prob",
          chrlabs = c(1:23, "Unplaced"), logp = T)


# Then a RDA ...


# RDA - ind ---------------------------------------------------------------

#library(codep)
#library(adespatial)
#library(adegraphics)
library(vegan)
library(ape)
library(car)

ID.Sag <- pop.data %>% dplyr::filter(Region_echantillonnage == "Saguenay") %>% pull(ID_GQ)

dist.gl <- dist(tab(gl.data[indNames(gl.data) %nin% ID.Sag,  ]), method = "euclidean", diag = T, upper = T)

#dist.all.gl   <- dist(tab(gl.all), method = "euclidean", diag = T, upper = T)

as.matrix(dist.gl)[1:5, 1:5] 
as.matrix(dist.gl) %>% nrow()


# Sort env variable as the distance matrix

env.data <- data.frame(ID_GQ = attr(dist.gl, "Labels")) %>% left_join(pop.data) %>% left_join(Env.data)


#db-Rda Response variables
##Pcoa on genetic distance matrix. Genetic distances are then in the multivariate space format, adapted to db-RDA
Pcoa <- pcoa(dist.gl)
Pcoa

#"There were no negative eigenvalues. No correction was applied"   
#If negative or null eigenvalues are produced, they need to be expluded
#The first axe explains 49% of the variation. We will keep them all to analyse all genetic variation

#Extract Pcoa principal components, which will be the response variable in the db-RDA
X <-Pcoa$vectors
#Look at genotypes distribution in relation to the first 2 Pcoa axes
plot(X[,1], X[,2])
#Clusters are clearly appearent

#explanatory variables
#Create a matrix with all expanatory variables, 16 MEMs, depth, and n-1 years coded in dummy variables (2001, 2002, 2008, 2013, 2014)  
Y <- env.data %>% dplyr::select(BO22_tempmean_bdmean, BO22_tempmean_ss, BO22_salinitymean_bdmean, BO22_salinitymean_ss, BO22_dissoxmean_bdmean, BO_bathymean)

plot(rda(Y))

#Correlation among explanatory variables
##Look at correlation among explanatory variables
cor(Y)
#Nothing problematic. Not surprising cause MEM are orthogonal to each other

#db-RDA global model with all explanatory variables
rda1 <- rda(X, Y)

#Looking at VIf for multicolinearity within the model. Greater than around 10 is problematic
vegan::vif.cca(rda1)
#We don't need to exclude variables because of correlation

#Explained variance by the global db-RDA model. Adjusted R2 accounts for the number of variables 
RsquareAdj(rda1)
#36% explained variance

#db-RDA Global model probability
anova(rda1, perm=99)
#p=0.001 

screeplot(rda1)
plot(rda1, scaling = 3)

#Variables selection with OrdiR2Step. 
#This function allows to add and remove variables to maximise the explained variance. 
#To avoid overfitting, selected variables should not explain more than the global model (36%) 

#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0<-rda(X ~ 1, Y)
#OrdiR2step will move towards the global model with all explanatory variables
rdaG<- rda(X ~ ., Y)

#Variables selection

#Variables will be selected until the 36% (rda global model) is reached
Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="both")                
#summary table with selected variables               
Sel$anova



Ysel <-  env.fas %>% dplyr::select(BO22_tempmean_ss, BO22_salinitymean_bdmean, BO22_tempmean_bdmean, Longueur_mm)
rdaS <- rda(X ,Ysel)

rdaS <- rda1

summary(rdaS, scaling=1)  


site.df <- data.frame(ID_GQ = dimnames(rdaS$Ybar)[[1]],
                       score = scores(rdaS, display="sites", choices=c(1,2), scaling=1)) %>%
  left_join(pop.data) 


arrow.df <- data.frame(name = scores(rdaS, display="bp", choices=1, scaling=1) %>% row.names(),
                       x0 = 0,
                       y0 = 0,
                       xmax = scores(rdaS, display="bp", choices=1, scaling=1) %>% as.vector(),
                       ymax = scores(rdaS, display="bp", choices=2, scaling=1) %>% as.vector())

graph.RDA <-  site.df %>% 
  left_join(ID_Ecoregion) %>% 
  ggplot(aes(x = score.RDA1, y = score.RDA2)) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  geom_point(aes(col = Ecoregion), cex = 2) +
  scale_colour_brewer(palette = "Set1") +
  #scale_colour_manual(values = c("orange", "khaki", "firebrick2", "royalblue2")) +
  #scale_fill_manual(values = c("ivory", "khaki","deeppink", "red", "firebrick2", "royalblue2")) +
  
  geom_segment(data = arrow.df , aes(x = x0, y = y0, xend = xmax*4, yend= ymax*4),
               arrow = arrow(), cex = 0.7)+
  geom_text(data = arrow.df , aes(label = str_remove(name, "BO22_"), x = xmax*2, y = ymax*2),
            nudge_x = + 0.05, nudge_y = -0.01) +
  
  labs(x=c("db-RDA-1"), y=c("db-RDA2")) + 
  annotate("text", label = paste("R²-adjusted =", RsquareAdj(rdaS)$adj.r.squared %>% as.numeric() %>%  round(digits = 5)),
           x = 0.3, y = 0.2, vjust = "inward", hjust = "inward") +
  theme_bw() 
graph.RDA

ggsave(filename = here::here("02_Results/02_PopulationGenetics/09_RDA/RDA_FAS_BioOracle.png"), 
       plot =graph.RDA.fas, 
       width =7, height =5 , units = "in",
       dpi = 300)



library(terra)
library(ggspatial)
library(tidyterra)


admin <- terra::vect(rnaturalearth::ne_countries(scale = "medium", continent = "north america", returnclass	  = "sf"))

extent(admin)

plot(terra::crop(env.bioOracle, extent(admin))



Set1.col <- RColorBrewer::brewer.pal(5,"Set1")


 Q.pie <- k3.df %>% left_join(pop.data) %>% 
  group_by(Numero_unique_groupe, Clust) %>%
  summarise(N = n(),
            Lat = mean(Latitude_echantillonnage_DD),
           Long = mean(Longitude_echantillonnage_DD)) %>% 
   pivot_wider(names_from = Clust, values_from = N, values_fill = 0, names_prefix = "c") %>% 
   mutate(Ntotal = c1+c2+c3)

 Q.pie

 clust.terra <-terra::vect(k3.df  %>% left_join(pop.data) %>%  
                             group_by(Numero_unique_groupe, Longitude_echantillonnage_DD, Latitude_echantillonnage_DD, Clust) %>% 
                             #    group_by(Numero_unique_groupe, Longitude_echantillonnage_DD, Latitude_echantillonnage_DD, Cat) %>%
                             summarize( value = n()) %>% mutate(N = sum(value)), geom = c("Longitude_echantillonnage_DD", "Latitude_echantillonnage_DD"), keepgeom = T, crs = "+proj=longlat")
 pop.terra <-terra::vect(pop.data %>%  
                             group_by(Numero_unique_groupe, Longitude_echantillonnage_DD, Latitude_echantillonnage_DD) %>% 
                             #    group_by(Numero_unique_groupe, Longitude_echantillonnage_DD, Latitude_echantillonnage_DD, Cat) %>%
                             summarize( value = n()) %>% mutate(N = sum(value)), geom = c("Longitude_echantillonnage_DD", "Latitude_echantillonnage_DD"), keepgeom = T, crs = "+proj=longlat")
 
 
 library(scatterpie)
gg.map <- ggplot() +
  geom_sf(data = admin, fill=NA, size=0.1) + 
  geom_sf(data = eco.shp[eco.shp$LABEL %in% c(6,8,9,10,12),], fill=NA, aes(col = NAME_E)) +
  geom_sf(data = clust.terra,
                    aes(cex = value)) +          
  #geom_scatterpie(aes(x=Long, y=Lat, group =Numero_unique_groupe, r = Ntotal*20000),
  #                data = Q.pie, cols = c("c1","c2", "c3"),  long_format = F) +
  facet_grid(~Clust) +
  coord_sf( crs = sf::st_crs("EPSG:3573"),xlim = c(-2000000,3000000), ylim = c(-4500000, -1000000)) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11),
        legend.position = "none")
gg.map

plot(env.bioOracle) +   coord_sf( crs = sf::st_crs("EPSG:3573"),xlim = c(1000000,3000000), ylim = c(-3500000, -2000000)) 
  
plot(pop.terra, add= T)

library(rasterVis)
library(raster)

plot(crop(env.bioOracle, extent(-130, -55, 45, 75)))

testy <- terra::rast(crop(env.bioOracle, extent(-75, -50, 52, 65)))

plot(crop(env.bioOracle, extent(-78, -50, 55, 65)))


plot(testy$BO22_tempmean_bdmean)
plot(pop.terra, add= T)

plot(testy$BO_bathymean)
plot(pop.terra, add= T)
     

plot(testy$BO22_tempmean_ss)
plot(pop.terra, add= T)


test<- terra::rast(env.bioOracle) 

plot(test)

gg.map <- ggplot() +
 # geom_raster(data = test, aes(fill = BO22_tempmean_bdmean))+
  geom_sf(data = admin, fill=NA, size=0.1) + 
  geom_sf(data = eco.shp[eco.shp$LABEL %in% c(6,8,9,10,12),], fill=NA, aes(col = NAME_E)) +
  geom_sf(data = clust.terra,
          aes(cex = value)) +          
  #geom_scatterpie(aes(x=Long, y=Lat, group =Numero_unique_groupe, r = Ntotal*20000),
  #                data = Q.pie, cols = c("c1","c2", "c3"),  long_format = F) +
  facet_grid(~Clust) +
  coord_sf( crs = sf::st_crs("EPSG:3573"),xlim = c(1000000,3000000), ylim = c(-3500000, -2000000)) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

gg.map







