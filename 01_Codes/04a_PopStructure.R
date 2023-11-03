# Info --------------------------------------------------------------------

# Population genetics analysis of Boreogadus saida
#  at global spatial scale
#  540 individuals
#
# Audrey Bourret
# 2023-08-22
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




pop.data %>% mutate(date_extraction = paste(Jour_extraction, Mois_extraction, Annee_extraction, sep = "-")) %>% 
             group_by(Numero_unique_groupe, date_extraction, Responsable_extraction) %>% 
             summarise(Ntotal = n()) %>% 
             dplyr::filter(Responsable_extraction == "S_Khan") %>% 
             View()
 

# Genetic Data ------------------------------------------------------------

load( file.path("./00_Data/06b_Filtering.ref", "A_ALL_samples", "07_Final", "populations.28309snps.540ind.adegenet.Rdata"))


gl.final
gi.final


na.gi.count <- function(gi){
  res <- apply(tab(gi), MARGIN = 2, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  res <- res[str_ends(names(res), "[.]0")] 
  
  names(res) <- names(res) %>% str_remove("[.]0")
  
  return(res)
  
}

# Function to create a list of loci, from a genind object

filter.MAF.NA <- function(gi, MAF.trs = 0.5, NA.trs = 0.5){
  # Create vectors for each loci
  MAF.res <- adegenet::minorAllele(gi)
  NA.res  <- na.gi.count(gi)
  
  # Filter by threshold
  MAF.loc <- dimnames(MAF.res[MAF.res >= MAF.trs])[[1]]
  cat("There is", length( MAF.loc), "loci with MAF =", MAF.trs, "\n")
  
  NA.loc <- names(NA.res[NA.res <= NA.trs])
  cat("There is", length(NA.loc), "loci with NA =", NA.trs, "\n")
  
  # LOCI with both conditions
  LOCI.res <- c(MAF.loc, NA.loc)[duplicated(c(MAF.loc, NA.loc)) == T]
  LOCI.res %>% length()
  
  cat("There is", length(LOCI.res), "loci with BOTH MAF =", MAF.trs, "and NA =" , NA.trs, "\n")
  
  return(LOCI.res)
}

LOC.MAF05.NA10 <- filter.MAF.NA(gi.final, MAF.trs = 0.05, NA.trs = 0.10)
LOC.MAF10.NA05 <- filter.MAF.NA(gi.final, MAF.trs = 0.10, NA.trs = 0.05)



# PCA ---------------------------------------------------------------------

pca.MAF10.NA05  <- glPca(gl.final[, locNames(gl.final) %in% LOC.MAF10.NA05], center = TRUE, scale = FALSE,  
                   parallel = TRUE, n.core =8, nf = 1000)


#save(list = c("pca.test" ),
#     file = here("02_Results/01_PopStruct", "01_PCA", "PCA.Rdata"))

load( here::here("02_Results/01_PopStruct", "01_PCA", "PCA.Rdata"))

QuickPop::pca_varplot(pca.test)

gPCA.extraction1 <- pca.MAF10.NA05 %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Responsable_extraction)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~Region_echantillonnage, ncol = 3) +
  #stat_ellipse(aes(col = Numero_unique_groupe))+
  geom_point(alpha = 0.5, size = 2) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) #+ theme(legend.position = "none")
gPCA.extraction1



gPCA.extraction2 <- pca.MAF10.NA05 %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  dplyr::filter(Responsable_extraction == "S_Khan") %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Region_echantillonnage)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~paste(Mois_extraction, Jour_extraction), ncol = 3) +
  #stat_ellipse(aes(col = Numero_unique_groupe))+
  geom_point(alpha = 0.5, size = 2) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8)
gPCA.extraction2


gPCA.extraction <- ggpubr::ggarrange(gPCA.extraction1 + ggtitle("Reponsables") , 
                                     gPCA.extraction2 + ggtitle("S_Khan, by day"))
gPCA.extraction

ggsave(filename = file.path(here::here(), "02_Results", "01_PopStruct", "01_PCA", "PCA_Extraction.png"), 
       plot = gPCA.extraction,
       height = 6, width = 10, units = "in", bg = "white")   



gPCA.PC12 <- pca.MAF10.NA05 %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate(CAT = ifelse(Region_echantillonnage %in% c("0A", "0B", "1A"), "01AB",
               ifelse(Region_echantillonnage %in% c("2G", "2H", "2J"), "2GHJ",
               ifelse(Region_echantillonnage %in% c("3K", "4S"), "3K4S",        
                      Region_echantillonnage)))) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = CAT)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(~Region_echantillonnage, ncol = 3) +
  stat_ellipse()+
  geom_point(alpha = 0.5, size = 2) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) #+ theme(legend.position = "none")
gPCA.PC12

gPCA.PC34 <- pca.MAF10.NA05 %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate(CAT = ifelse(Region_echantillonnage %in% c("0A", "0B", "1A"), "01AB",
                      ifelse(Region_echantillonnage %in% c("2G", "2H", "2J"), "2GHJ",
                             ifelse(Region_echantillonnage %in% c("3K", "4S"), "3K4S",        
                                    Region_echantillonnage)))) %>% 
  ggplot(aes(x = score.PC3, y = score.PC4, col = CAT)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(~Region_echantillonnage, ncol = 3) +
  stat_ellipse()+
  geom_point(alpha = 0.5, size = 2) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC3 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[3] %>% round(3) *100, "%)"),
    y = paste0("PC4 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[4] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) #+ theme(legend.position = "none")
gPCA.PC34


gPCA.PC56 <- pca.MAF10.NA05 %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate(CAT = ifelse(Region_echantillonnage %in% c("0A", "0B", "1A"), "01AB",
                      ifelse(Region_echantillonnage %in% c("2G", "2H", "2J"), "2GHJ",
                             ifelse(Region_echantillonnage %in% c("3K", "4S"), "3K4S",        
                                    Region_echantillonnage)))) %>% 
  ggplot(aes(x = score.PC5, y = score.PC6, col = CAT)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(~Region_echantillonnage, ncol = 3) +
  stat_ellipse()+
  geom_point(alpha = 0.5, size = 2) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC5 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[5] %>% round(3) *100, "%)"),
    y = paste0("PC6 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[6] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) #+ theme(legend.position = "none")
gPCA.PC56



gPCA <- ggpubr::ggarrange(gPCA.PC12 , 
                         gPCA.PC34, common.legend = T)
gPCA

ggsave(filename = file.path(here::here(), "02_Results", "01_PopStruct", "01_PCA", "PCA.png"), 
       plot = gPCA,
       height = 6, width = 10, units = "in", bg = "white")   



pop.data %>% View()

pca.MAF10.NA05 %>% QuickPop::pca_scoretable(naxe = 6) %>%
       left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
   dplyr::filter(Region_echantillonnage == "Beaufort_Sea_ecoregion", score.PC2 > -2) %>% View()





count.ind.na.gl <- function(gl){
  res <- apply(tab(gl,  NA.method = c("asis")), MARGIN = 1, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  return(res)
  
}


na.info <- data.frame(ID_GQ = indNames(gl.final),
                      NNA = count.ind.na.gl(gl.final[, locNames(gl.final) %in% LOC.MAF10.NA05]))




na.info %>% left_join(pop.data) %>% 
  ggplot(aes(x = str_remove(Region_echantillonnage, "_ecoregion"), y = NNA, fill = NNA)) +
  geom_boxplot() + geom_jitter(pch = 21) +
  scale_fill_distiller(palette = "Spectral") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

gPCA.NA <- pca.MAF10.NA05 %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(na.info, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = NNA)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~Region_echantillonnage) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  scale_color_distiller(palette = "Spectral") +
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.MAF10.NA05)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA.NA


ggsave(filename = file.path(here::here(), "02_Results", "01_PopStruct", "01_PCA", "PCA_wNA.png"), 
       plot = gPCA.NA,
       height = 6, width = 6, units = "in", bg = "white")   




# Admixture ---------------------------------------------------------------


current.wd <- getwd()

# Create a BED files ...

# CREATE VCF WITH UNIQUE

vcf.path <- file.path(here::here(), "00_Data/06b_Filtering.ref/A_ALL_samples/07_Final/populations.28309snps.540ind.H06.DP.single.final.recode.vcf")

data.frame(ID = LOC.MAF05.NA10)

write.csv(data.frame(ID = LOC.MAF05.NA10), file.path( "02_Results/01_PopStruct/02_Admixture", "Loc.MAF05.NA10.csv"),
          row.names = F, quote = F)

cmd <- paste("--vcf", vcf.path,
             "--recode",
             #paste("--indv",ID.fas, collapse = " "),
             "--snps",file.path("02_Results/01_PopStruct/02_Admixture", "Loc.MAF05.NA10.csv"),
             "--out", file.path("02_Results/01_PopStruct/02_Admixture", paste0("Bsaida.", length(LOC.MAF05.NA10),"snps.540ind.final"))
)

cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)

tail(A)

# A2

cat(file =  file.path(current.wd, "02_Results/01_PopStruct/02_Admixture", "VCF.fas.filter.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")

cmd1 <- paste("--vcf", file.path(current.wd, "02_Results/01_PopStruct/02_Admixture", "Bsaida.18531snps.540ind.final.recode.vcf" ),
              #"--recode",
              "--plink-tped",
              "--out", file.path(current.wd, "02_Results/01_PopStruct/02_Admixture", "Bsaida.18531.540ind.final"))


cmd1

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)
A1

cmd2a <- paste("--tfam", "./02_Results/01_PopStruct/02_Admixture/Bsaida.18531.540ind.final.tfam",
               "--tped", "./02_Results/01_PopStruct/02_Admixture/Bsaida.18531.540ind.final.tped",
               "--make-bed",
               "--out", "./02_Results/01_PopStruct/02_Admixture/Bsaida.18531.540ind.final")



A2a <- system2("/home/genyoda/Documents/Programs/plink_linux_x86_64_20210606/plink", cmd2a, stdout=T, stderr=T)
A2a

# Original SNP

bed.file <- file.path(here::here(), "./02_Results/01_PopStruct/02_Admixture/Bsaida.18531.540ind.final.bed" )
file.exists(bed.file)
fam.file <- bed.file %>% str_replace(".bed", ".fam")
fam <- read.table(fam.file)


for(k in 1:10){
  
  print(k)  
  
  setwd(file.path(here::here(), "/02_Results/01_PopStruct/02_Admixture/") ) 
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               bed.file,
               k, # the number of K
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("k",k, ".log"),
      "\n", cmd, "\n",
      A, # what to put in my file
      append= F, sep = "\n")
  
  setwd(here::here())
  
}

# Cross-validation results:

CV.res <- data.frame(k = 1:10,
                         CV = NA,
                         stringsAsFactors = F)


for(i in 1:nrow(CV.res)){
  # Which k
  k <- CV.res[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path("./02_Results/01_PopStruct/02_Admixture/", paste0("k",k, ".log")))
  CV.temp <- temp %>% str_subset("CV error")
  CV <- sapply(str_split(CV.temp, ":"), `[`, 2) %>% str_remove_all(" ")
  
  # Add to my data.frame
  CV.res[i, "CV"] <- CV
  
}

CV.res$CV <- as.numeric(as.character(CV.res$CV))

CV.res %>% arrange(CV)

plot(CV.res$CV)

gg.CV <- CV.res %>% mutate(color = ifelse(k == 4, "red", "black")) %>% 
  ggplot(aes(x = factor(k), y = CV)) + 
  geom_point(size = 2, aes(col = color)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "K", y = "Cross-validation error") +
  theme_bw() +
  theme(legend.position = "none")
gg.CV

ggsave(filename = file.path(here::here(), "02_Results/01_PopStruct/02_Admixture", "Admixture.CV.png"), 
       plot = gg.fas,
       height = 3.5, width = 4, units = "in")   


k <- 5

Q.k2.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/02_Admixture", paste0("Bsaida.18531.540ind.final.",2,".Q")))
Q.k3.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/02_Admixture", paste0("Bsaida.18531.540ind.final.",3,".Q")))
Q.k4.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/02_Admixture", paste0("Bsaida.18531.540ind.final.",4,".Q")))
Q.k5.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/02_Admixture", paste0("Bsaida.18531.540ind.final.",5,".Q")))

Q.res <- bind_rows(cbind(fam$V1, Q.k5.res, K = 5),
                   cbind(fam$V1, Q.k4.res, K = 4),
                   cbind(fam$V1, Q.k3.res, K = 3),
                   cbind(fam$V1, Q.k2.res, K = 2))


head(Q.res)

#Q.fas.res <- cbind(fam.fas$V1, Q.fas.res)

names(Q.res) <- c("ID_GQ", paste0("Q", 1:k), "K")

#reorder(ID, Qvalue, FUN = function(x) min(x))

gg.str <- Q.res %>% pivot_longer(cols =  paste0("Q", 1:k), names_to = "Group", values_to = "Q") %>% 
  mutate(Group = factor(Group, levels = c("Q2", "Q4", "Q3", "Q5", "Q1"))) %>% 
  left_join(pop.data) %>% 
  #left_join(clust.k5.df ) %>%
  #left_join(na.info) %>% 
  #mutate(CAT = ifelse(NNA > 0.05, "High", "Low")) %>% 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  
  geom_col() +
  #facet_grid(. ~Lieu_echantillonnage + Mois_echantillonnage, space = "free", scale = "free") +
  facet_grid(K ~ CAT + Clust , space = "free", scale = "free") +
  scale_fill_brewer(palette = "Set1") +
  labs(y="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )
gg.str




clust.terra <-terra::vect(clust.k5.df  %>% left_join(pop.data) %>%  
                             group_by(Numero_unique_groupe, Longitude_echantillonnage_DD, Latitude_echantillonnage_DD, Clust) %>% 
                           #    group_by(Numero_unique_groupe, Longitude_echantillonnage_DD, Latitude_echantillonnage_DD, Cat) %>%
                           summarize( value = n()) %>% mutate(N = sum(value)), geom = c("Longitude_echantillonnage_DD", "Latitude_echantillonnage_DD"), keepgeom = T, crs = "+proj=longlat")



library(terra)
library(ggspatial)
library(tidyterra)


admin <- terra::vect(rnaturalearth::ne_countries(scale = "medium", continent = "north america", returnclass	  = "sf"))

Set1.col <- RColorBrewer::brewer.pal(5,"Set1")

gg.str.map <- ggplot(clust.terra) + 
  geom_sf(data = admin, fill=NA, size=0.1) +
  geom_sf(data = clust.terra,
          aes(col = Clust, cex = value)) +
   #coord_sf( crs = sf::st_crs("EPSG:6622"),xlim = c(-5000000,1000000), ylim = c(-50000, 4500000)) +
  coord_sf( crs = sf::st_crs("EPSG:3573"),xlim = c(-2000000,3000000), ylim = c(-4500000, -1000000)) +
    #coord_sf(xlim = c(-140, -54), ylim = c(40, 80), expand = FALSE) +
  facet_grid(Clust ~ .) +
  #facet_wrap(.~Clust) +
  scale_color_manual(values = c(Set1.col[c(1,3,4)], "brown", Set1.col[5])) +
  #scale_colour_brewer(palette = "Set1") +     
  xlab("Longitude") + ylab("Latitude") +
  theme_bw()
gg.str.map


gg.str.complete <-  ggarrange(gg.str + theme(legend.position = "none"), 
          gg.str.map,
          widths = c(3,2))


ggsave(filename = file.path(here::here(), "02_Results", "01_PopStruct", "Structure_wMap.png"), 
                  plot = gg.str.complete,
                  height = 6, width = 10, units = "in", bg = "white")



# SnapClust ---------------------------------------------------------------

library(adegenet)


test.clust <- snapclust.choose.k(x = gi.final[loc = LOC.MAF05.NA10],
                                 max = 10)


plot(test.clust, type = "b", cex = 2, xlab = "k", ylab = "AIC")
points(which.min(test.clust), min(test.clust), col = "blue", pch = 20, cex = 2)
abline(v = 5, lty = 2, col = "red")

clust.k5 <- snapclust( gi.final[loc = LOC.MAF05.NA10], k = 5)

clust.k5.prob <-  data.frame(ID_GQ = dimnames(clust.k5$proba)[[1]], 
                                 Q1 = clust.k5$proba[,1],
                                 Q2 = clust.k5$proba[,2],
                                 Q3 = clust.k5$proba[,3],
                                 Q4 = clust.k5$proba[,4],
                                 Q5 = clust.k5$proba[,5]) %>% 
  mutate(Q1.corr = ifelse(is.na(Q1), 1 - (Q2 +Q3+ Q4 + Q5), Q1),
         Q2.corr = ifelse(is.na(Q2), 1 - (Q1 +Q3+ Q4 + Q5), Q2),
         Q3.corr = ifelse(is.na(Q3), 1 - (Q1 +Q2+ Q4 + Q5), Q3),
         Q4.corr = ifelse(is.na(Q4), 1 - (Q1 +Q2+ Q3 + Q5), Q4),
         Q5.corr = ifelse(is.na(Q5), 1 - (Q1 +Q2+ Q3 + Q4), Q5)
  ) %>% 
  dplyr::select(ID_GQ, Q1 = Q1.corr, Q2 = Q2.corr, Q3 = Q3.corr, Q4 = Q4.corr, Q5 = Q5.corr) %>% 
  pivot_longer(-ID_GQ, names_to = "Group", values_to = "Q") %>% 
  mutate(Method = "SnapClust")

clust.k5$proba %>% min(na.rm = T)

compoplot(clust.k5)

#save(list = c("test.clust", "clust.k5", "clust.k5.prob"),
#     file = "02_Results/01_PopStruct/03_SnapClust/FastClust.RData")

clust.k5.df <- data.frame(ID_GQ =  names(clust.k5$group),  
                          Clust = clust.k5$group)



# Fst ---------------------------------------------------------------------



table.fst <- function(fst){
  res <-  fst$Bootstraps %>% dplyr::select(Population1, Population2, "Lower bound CI limit", "Upper bound CI limit", "p-value", "Fst")
  
  return(res)
  
}

heat.fst <- function(fst){
  res <- bind_rows(table.fst(fst),
                   table.fst(fst) %>% mutate(Pop1.int = Population2,
                                             Pop2.int = Population1,
                                             Population1 = Pop1.int,
                                             Population2 = Pop2.int) %>% 
                     dplyr::select(-c(Pop1.int, Pop2.int))
  )
  
  return(res)         
  
}


pop.data %>% group_by(Annee_echantillonnage, Mois_echantillonnage, Region_echantillonnage) %>% 
  summarise(N = n()) %>% arrange(desc(N)) %>% 
  ggplot(aes(x = Mois_echantillonnage, y = N, fill = Region_echantillonnage)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Annee_echantillonnage)



pop(gl.final) <-  data.frame(ID_GQ = indNames(gl.final)) %>% 
  left_join(pop.data) %>% pull(Region_echantillonnage)

pop(test.gl) <-  data.frame(ID_GQ = indNames(test.gl)) %>% 
  left_join(ref.all.data) %>% pull(CAT)


table(pop(gl.final))



fst <- dartR::gl.fst.pop(gl.final[,locNames(gl.final) %in% LOC.MAF05.NA10],
                         nboots = 99, percent = 95, nclusters = 20)
ID.test <- na.info %>% dplyr::filter(NNA <=.1) %>% pull(ID_GQ)
fst.NAind10 <- dartR::gl.fst.pop(gl.final[indNames(gl.final) %in%  ID.test,locNames(gl.final) %in% LOC.MAF05.NA10],
                         nboots = 99, percent = 95, nclusters = 20)



names(fst)

fst$Bootstraps %>% dplyr::select(Population1, Population2, `p-value`,Fst) %>% arrange(desc(Fst))

fst$Bootstraps$Fst %>% max()


fst %>% heat.fst() %>%  #  left_join(pop.data %>% dplyr::select(Population1 = Numero_unique_groupe, Region1 = Region_echantillonnage, Mois1 = Mois_echantillonnage) %>% distinct()  ) %>% 
  #left_join(pop.data %>% dplyr::select(Population2 = Numero_unique_groupe, Region2 = Region_echantillonnage, Mois2 = Mois_echantillonnage) %>% distinct()  ) %>% 
  #left_join(pop.data %>% dplyr::select(Population2 = RegionAssesment, SFA2 = RegionAssesment)  %>% distinct()) %>% #pull(SFA1) %>% unique()
  mutate(   Sign = ifelse(`p-value` <= 0.05, "*", ""),
         Population1 = str_remove(Population1, "_ecoregion"),
         Population2 = str_remove(Population2, "_ecoregion"),
         Population1 = factor(Population1, levels = c("Beaufort_Sea", "Eastern_Arctic",  "1A", "0A", "0B", "WAZ", "Hudson_Bay", "2G", "2H", "2J", "3K", "4S", "Saguenay")),
                              
         Population2 = factor(Population2, levels = c("Beaufort_Sea", "Eastern_Arctic",  "1A", "0A", "0B", "WAZ", "Hudson_Bay", "2G", "2H", "2J", "3K", "4S", "Saguenay"))
         
           ) %>% 
  #dplyr::filter(Mois1 == "9", Mois2 == "9") %>%              
  #  dplyr::filter(Region1 == "4X", Region2 == "4X") %>%      
  ggplot(aes(x=Population1, y=Population2, fill=Fst)) + 
  geom_tile(colour=  "gray") +
  # geom_point(aes(x = Population1, y=Population2), size = 3)+
  geom_text(aes(label = Sign), cex = 3) +
  #scale_y_discrete(limits=rev) +
  #scale_x_discrete(limits=rev) +
  scale_fill_gradient(low = "ivory1", high = "red3", na.value = "white", limits = c(0,0.05))+
  # scale_fill_gradient(low = "ivory1", high = "dodgerblue2")+
  scale_shape_manual(values = c("","*"), guide = "none") +
  #scale_fill_distiller(palette = "Spectral") +
  labs(x = NULL, y = NULL) +
  #theme_bw() +
 # facet_grid(Region1 ~ Region2, space = "free", scale = "free") +
  theme_minimal() + 
  theme(#axis.text.x = element_blank(),
    strip.text = element_text(angle = 0),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "cm"),
    panel.border = element_rect(fill = NA, colour = "black"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    #legend.background = element_rect(colour = "black", fill = "grey95", size = 0.2),
    #legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle =  90, vjust = 0.5, hjust = 1))



# He vs H0 ----------------------------------------------------------------

pop(gi.final) <-  data.frame(ID_GQ = indNames(gi.final)) %>% 
  left_join(pop.data) %>% pull(Region_echantillonnage)

basic_stat_Region <- hierfstat::basic.stats(gi.final[loc = LOC.MAF05.NA10], diploid = TRUE)


gi.test <- gi.final[loc = LOC.MAF05.NA10]
gi.test <- gi.test[ID.test,]
basic_stat_Region.test <- hierfstat::basic.stats(gi.test, diploid = TRUE)

#save(list = c("basic_stat_Region"), 
#     file = file.path(here::here(), "02_Results/01_PopStruct/05_HeHo/HeHo.RData"))

load(file.path(here::here(), "02_Results/01_PopStruct/05_HeHo/HeHo.RData"))

# Homemade function to extract important infos
extract.basic.stats <- function(basic_stat_ALL){
  data.frame(ID = dimnames(basic_stat_ALL$Fis)[[2]],
             Fis = apply(basic_stat_ALL$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE),
             Fis.IC95.low = apply(basic_stat_ALL$Fis, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[1]),
             Fis.IC95.high = apply(basic_stat_ALL$Fis, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[2]),
             
             He = apply(basic_stat_ALL$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE),
             He.IC95.low = apply(basic_stat_ALL$Hs, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[1]),
             He.IC95.high = apply(basic_stat_ALL$Hs, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[2]),
             
             Ho = apply(basic_stat_ALL$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE),
             Ho.IC95.low = apply(basic_stat_ALL$Ho, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[1]),
             Ho.IC95.high = apply(basic_stat_ALL$Ho, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[2])
             
  )
}

Hstat_Region      <- extract.basic.stats(basic_stat_Region)
#Hstat_Region      <- extract.basic.stats(basic_stat_Region.test)

#write_csv(Hstat_Region , file.path(here::here(), "02_Results/01_Overall_PopGen/00_BasicStats/Hstat_Region.csv"))
#write_csv(Hstat_Gen_ZONE , file.path(here::here(), "02_Results/01_Overall_PopGen/00_BasicStats/Hstat_Gen_ZONE.csv"))
#write_csv(Hstat_Gen_ZONE_FG , file.path(here::here(), "02_Results/01_Overall_PopGen/00_BasicStats/Hstat_Gen_ZONE_FG.csv"))


Hstat_Region

het.gg.region <- rbind(Hstat_Region %>% dplyr::select(ID, value = Ho, IC95.high = Ho.IC95.high, IC95.low =  Ho.IC95.low ) %>%
                         mutate(name = "Ho"),
                       Hstat_Region %>% dplyr::select(ID, value = He, IC95.high = He.IC95.high, IC95.low =  He.IC95.low ) %>%
                         mutate(name = "He") ) %>% 
  
  ggplot(aes(x = ID, y = value, fill = name, group = name, shape = name)) + 
  #scale_y_continuous(limits = c(0.15, 0.20)) +
  
  geom_errorbar(aes(ymin=IC95.low, ymax=IC95.high), width=0, position=position_dodge(width=0.3)) +
  geom_point(position=position_dodge(width=0.3)) + 
  scale_fill_manual(values = c("black", "darkgray")) +
  scale_shape_manual(values = c(24,25)) +
  labs(y= "Heterozygosity") +
  #facet_grid(~RegionAssesment, space = "free", scale = "free") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
het.gg.region
