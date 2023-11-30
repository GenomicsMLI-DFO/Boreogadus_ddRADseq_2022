# Info --------------------------------------------------------------------

# PCADAPT for Boreogadus saida
#  at global spatial scale
#  549 individuals (9 Arctogadus, 45 intermediate, 500 Boreogadus)
#
# Audrey Bourret
# 2023-10-27
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

vcf.path <- "/media/genyoda/Extra_Storage/Projets/Data_Trevor/02_Results/08_Boreogadus/10x/populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.noinver.recode.vcf"

vcf.data <- vcfR::read.vcfR(vcf.path)
gl.data  <- vcfR::vcfR2genlight(vcf.data) 


cmd1 <- paste("--vcf", vcf.path, 
              #"--recode",
              "--plink-tped",
              "--out",  "02_Results/02_PCADAPT/populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.noinver")


cmd1

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)

cmd2a <- paste("--tfam", "./02_Results/02_PCADAPT/populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.noinver.tfam", 
               "--tped", "./02_Results/02_PCADAPT/populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.noinver.tped", 
               "--make-bed", 
               "--out", "./02_Results/02_PCADAPT/populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.noinver" 
               
)

 A2a <- system2("/home/genyoda/Documents/Programs/plink_linux_x86_64_20210606/plink", cmd2a, stdout=T, stderr=T)
A2a


# PCADAPT -----------------------------------------------------------------

library(pcadapt)
#BiocManager::install("qvalue")
library("qvalue")

# Convertion to plink .bed format

# Read .bed in PCAadapt
pcadapt.genotype  <- read.pcadapt(file.path("./02_Results/02_PCADAPT/populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.noinver.bed" ),
                                         type = "bed")


pcadapt.snp <- read.delim(file.path(here::here(), "./02_Results/02_PCADAPT/populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.noinver.bim" ),
                                 header = F) %>% pull(V2)

# Run pcadapt

K.init <- 10

pcadapt.k10   <- pcadapt(pcadapt.genotype, K =K.init)
# Check screeplot

plot(pcadapt.k10, option = "screeplot") 

# Check structure

plot(pcadapt.k10, option = "scores") 

# K = 2 pour sfa et K = 3 pour final

pcadapt.k2 <- pcadapt(pcadapt.genotype , K = 2)
pcadapt.k3 <- pcadapt(pcadapt.genotype , K = 3)
pcadapt.k4 <- pcadapt(pcadapt.genotype , K = 4)


plot(pcadapt.k2, option = "manhattan")

hist(pcadapt.k2$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

plot(pcadapt.k2, option = "qqplot")

plot(pcadapt.k2, option = "stat.distribution")

# Statistics
#x$pvalues 
alpha <- 0.05 

qval.k2 <- qvalue::qvalue(pcadapt.k2$pvalues)$qvalues
outliers.k2 <- which(qval.k2 < alpha)
length(outliers.k2)


padj.k2 <- p.adjust(pcadapt.k2$pvalues,method="BH")
outliers.padj.k2 <- which(padj.k2 < alpha)
length(outliers.padj.k2)

bonf.k2 <- p.adjust(pcadapt.k2$pvalues,method="bonferroni")
outliers.bonf.k2 <- which(bonf.k2 < alpha)
length(outliers.bonf.k2)


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

res.pcadapt %>% mutate(test = -log10(pvalues)) %>% arrange(pvalues) %>% View(
  
)


write_csv(res.pcadapt, "02_Results/02_PCADAPT/PCADAPT_k2_results.csv")

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

nLoc(gl.final)

library(RColorBrewer)
library(eulerr)

venn.pcadapt <- plot(euler(list(pcadapt.k2=   pcadapt.snp[outliers.k2]  , pcadapt.k3 =  pcadapt.snp[outliers.k3], pcadapt.k4 =  pcadapt.snp[outliers.k4] ) , shape = "circle"),
                     quantities = T,
                     fill = brewer.pal(n = 3, name = "Dark2"),
                     alpha = 0.7, edges = F)

venn.pcadapt




pcadapt.outliers <- pcadapt.snp[outliers.k2]



# Create new VCF ----------------------------------------------------------

pcadapt.outliers %>% length() / nLoc(gl.data)


SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, ","), `[`,1) %>% str_remove("scaffold"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )

SCAFFOLD.info

df.outlier <-  SCAFFOLD.info %>% 
  dplyr::filter(ID %in% pcadapt.outliers)

df.neutral <-  SCAFFOLD.info %>% 
  dplyr::filter(ID %nin% pcadapt.outliers)

nrow(df.outlier)
nrow(df.neutral)


write.csv(df.outlier %>% select(ID), file.path("02_Results/02_PCADAPT/Outlier_pcadapt_k2.csv"), 
          row.names = F, quote = F)

write.csv(df.neutral %>% select(ID), file.path("02_Results/02_PCADAPT/Neutral_pcadapt_k2.csv"), 
          row.names = F, quote = F)

# CREATE VCF WITH UNIQUE

vcf.path 

cmd <- paste("--vcf", vcf.path, 
             "--recode",
             "--snps", file.path("02_Results/02_PCADAPT/Outlier_pcadapt_k2.csv"),
             
             "--out", file.path("02_Results/02_PCADAPT/populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.noinver.outlier")
             )


cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)


cmd <- paste("--vcf", vcf.path, 
             "--recode",
             "--snps", file.path("02_Results/02_PCADAPT/Neutral_pcadapt_k2.csv"),
             
             "--out", file.path("02_Results/02_PCADAPT/populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.noinver.neutral")
)


cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)



# PCA ---------------------------------------------------------------------

library(Hmisc)

pca.outliers  <- glPca(gl.data[,locNames(gl.data) %in%  pcadapt.snp[c(outliers.bonf.k2) %>% unique()]], center = TRUE, scale = FALSE,  
                  parallel = TRUE, n.core =16, nf = 1000)

pca.neutral  <- glPca(gl.data[,locNames(gl.data) %nin% pcadapt.outliers], center = TRUE, scale = FALSE,  
                       parallel = TRUE, n.core =16, nf = 1000)


pca.outliers %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate( Region = ifelse(Region_echantillonnage %in% c("0A", "0B", "1A"), "0AB1A",
                                  ifelse(Region_echantillonnage %in% c("2G", "2H", "2J"), "2GHJ", 
                                         Region_echantillonnage ) 
  )) %>% 
  ggplot(aes(x = score.PC3, y = score.PC4, col = Region)) +
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


pca.neutral %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate( Region = ifelse(Region_echantillonnage %in% c("0A", "0B", "1A"), "0AB1A",
                          ifelse(Region_echantillonnage %in% c("2G", "2H", "2J"), "2GHJ", 
                                 Region_echantillonnage ) 
  )) %>% 
  left_join(Env.data) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = BO22_tempmean_ss)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(~Region_echantillonnage, ncol = 3) +
  #stat_ellipse(aes(col = Region))+
  geom_point(alpha = 0.5, size = 2) +  
  scale_color_viridis_c()+
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.neutral)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.neutral)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8)


