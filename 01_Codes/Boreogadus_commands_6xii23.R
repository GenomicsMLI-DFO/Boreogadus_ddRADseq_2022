# Info --------------------------------------------------------------------
#
# Overview: Analysis of ddRAD datasets in Boreogadus saida (Arctic cod)
# 
# Author: Audrey Bourret and Trevor Bringloe
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory
# Location: Maurice Lamontagne Institute
# Date: 2023-11-20
#

# Library -----------------------------------------------------------------
library(terra)
library(ggspatial)
library(tidyterra)
library(adegenet)
library(vcfR)
library(ggplot2)
library(viridis)
library(rnaturalearth)

# Map of sampling locations -----------------------------------------------------------------
#Import map data
world <- ne_countries(scale = "medium", returnclass = "sf")
#Read in data, tab delimited file with Lat, Long, Region, and Species columns
#Note, in order to plot Arctogadus records as seperate colour, these were integrated into Regions column
boreogadus.coord <- read.csv("Boreogadus_coords_inversions_21xii23.txt", sep = "\t")
#Order data by region
boreogadus.coord$Ecoregion <- factor(boreogadus.coord$Ecoregion, levels = c("Western Arctic", "Eastern Arctic", "Hudson Bay Complex", "Newfoundland-Labrador Shelves","Gulf of Saint Lawrence", "Arctogadus", "Boreogadus(Arctogadus MtDNA)", "Potential hybrid"))
#Create data file with GPS coordinates
boreogadus.coord.map <- boreogadus.coord %>% 
  filter(Chrom9Inversion %in% c("Heterozygous")) %>%
  terra::vect(geom = c("LonM", "LatM"),  crs = "+proj=longlat")

#Plot coordinates on map, transforming data to polar projection
#Regions depicted in viridis
gg.env <- ggplot(boreogadus.coord.map) + 
  geom_sf(data = world, fill="gray95", size=0.5, alpha = 1) +
  geom_sf(data = boreogadus.coord.map, alpha = .8) + 
  #scale_size_continuous(limits=c(0,25),breaks=c(1, 5, 10, 15, 20, 25)) + 
  #coord_sf(xlim = c(-2750000, 1250000), ylim = c(350000, 3750000), crs = sf::st_crs("EPSG:6622")) +
  coord_sf(crs = "+proj=laea +lat_0=90 +lon_0=-100 +x_0=40 +y_0=0") +
  scale_colour_manual(values=c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725","darkgrey", "white", "red3")) +
  #scale_size_manual(values=c("6", "2", "2", "2")) +
  scale_shape_manual(values=c(18, 19,4, 3)) +
  #scale_size_manual(values=c(1, 2, 0.5, 0.5)) +
  # Map limits
  #scale_shape_manual(values = c(16,15,17), name = "Year") +
  #coord_sf(xlim = c(-750000, 1250000), ylim = c(350000, 3750000), crs = sf::st_crs("EPSG:6622")) +
  # Others
  #facet_grid(. ~ pred.pop) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Beluga lcWGS") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), strip.text = element_text(size=10),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom")  
gg.env

ggsave(filename = file.path("boreogadus_hetero_regions.pdf"), 
       plot = gg.env + ggtitle("boreogadus regions"), #+ guides(colour = FALSE),
       height = 15, width = 24, units = "in")

# PCA and DAPC analyses -----------------------------------------------------------------
#Ran code for species and population level datasets
#Import vcf
boreogadus.pop.vcf <- read.vcfR("populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.chrom9inver.recode.vcf")
boreogadus.pop.genind <- vcfR2genind(boreogadus.pop.vcf)
boreogadus.pop.genlight <- vcfR2genlight(boreogadus.pop.vcf)

#Find clusters
grp.b.outlier <- find.clusters(boreogadus.pop.genlight, max.n.clust=20)
#Describe clusters using dapc
dapc.grp <- dapc(boreogadus.pop.genlight, boreogadus.pop.genlight$pop)

#Visualize dapc1
scatter(dapc.grp, 1, 2, col=c("#31688e", "#35b779", "#440154", "#fde725"))

#Create graph of loading values
#Change axes with axis option
contrib <- loadingplot(dapc.grp$var.contr, axis=1,
                       thres=.07, lab.jitter=1)

#plot graphs
pdf("boreogadus_dapc.group2.loadings.10x.pdf", height = 5, width = 8) 
#Plug in R code to create figure
dev.off()

#Export variable loadings
dapc.dataset.loadings.grp = as.data.frame(dapc.grp$var.contr)
dapc.dataset.pcaloadings.grp = as.data.frame(dapc.grp$pca.loadings)
dapc.dataset.da.grp = as.data.frame(dapc.grp$ind.coord)
dapc.dataset.tab.grp = as.data.frame(dapc.grp$tab)
dapc.dataset.grp = as.data.frame(dapc.grp$grp)
write.csv(dapc.dataset.loadings.grp, "dapc.dataset.daloadings.outlier.b.grp.csv")
write.csv(dapc.dataset.tab.grp, "dapc.dataset.tab.grp.csv")
write.csv(dapc.dataset.grp, "dapc.dataset.grp.csv")

#Regional, chromosomal, position information linked to datasets external to R
#Import data
dapc.dataset.loadings.grp <- read.csv("dapc.dataset.loadings.grp.csv", sep = ",")

#Run PCA
#Convert missing data to genotype means
x.boreogadus.pop <- tab(boreogadus.pop.genind, freq=TRUE, NA.method="mean")

#Run PCA, reading genind object as a genlight object
pca.boroegadus.chrom9.pop <- glPca(as.genlight(x.boreogadus.pop), center=TRUE)

#Create new dataframes of PCs and eigen values
pca.dataset.pop = as.data.frame(pca.boroegadus.chrom9.pop$scores)
pca.dataset.pop.eigen = as.data.frame(pca.boroegadus.chrom9.pop$eig)
#To be used in PCA plots
pve <- data.frame(PC = 1:496, pve = pca.dataset.pop.eigen/sum(pca.dataset.pop.eigen*100))
colnames(pve)[2] = "pve"

#Export dataframes
write.csv(pca.dataset.pop, "pca.dataset.pop.chrom9.csv")
#Import data back into R after adding Region column
pca.dataset.pop <- read.csv("pca.dataset.pop.chrom9.csv", sep = ",")
#Optional add group/population data to genind/genlight object
boreogadus.pop.genlight@pop <- as.factor(pca.dataset.pca.b.outlier.pop$Group)

#Reorder regional/MtDNA information
pca.dataset.species$MtDNA <- factor(pca.dataset.species$MtDNA, levels = c("Boreogadus saida", "Arctogadus glacialis"))
pca.dataset.pop$Ecoregion <- factor(pca.dataset.pop$Ecoregion, levels = c("Western Arctic", "Eastern Arctic", "Hudson Bay Complex", "Newfoundland-Labrador Shelves","Gulf of Saint Lawrence"))

#Plot PCAs with eclipses # plot species with mtDNA for col
gg.pca.pop <- ggplot(pca.dataset.pop, aes(x = PC1, y = PC2, col = Het.Chrom9)) + geom_point() +
  #scale_colour_manual(values=c("#440154", "#46327e", "#365c8d", "#277f8e", "#1fa187", "#4ac16d", "#a0da39", "#fde725")) +
  #scale_colour_manual(values=c("#440154", "#fde725")) +
  #scale_colour_manual(values=c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")) +
  scale_colour_viridis(option="A", limits = c(-0.50, 0.1), oob = scales::squish) +
  #scale_colour_viridis_c(direction=-1) +
  #scale_colour_manual(values=c("#2ed7e8ff", "darkblue", "darkgrey"))+
  scale_shape_manual(values=c(15, 16, 17, 18)) +
  guides(colour = guide_legend(title = "Inbreeding", order = 1),
         label.position = "right",
         title.position = "top", title.hjust = 0,
         override.aes = list(size = 5),
         order = 1) +
  theme(legend.position = c(0.75,0.45),
        legend.spacing.y = unit(0.25, "cm"),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 16),
        legend.background = element_blank()) +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  #xlab(paste0("PC1 (", signif(pve$pve[1], 1), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 2), "%)")) + 
  #optional stat eclipse
  #stat_ellipse() +
  theme_bw()
gg.pca.pop
gg.pca.boreogadus.neutral.pop

#print plot
ggsave(filename = file.path("boreogadus_dapc.chrom9.pca_1+2.het2.pop.10x.pdf"), 
       plot = gg.pca.pop + ggtitle("boreogadus dapc pop 10x  SNPs"), #+ guides(colour = FALSE) + guides(shape = FALSE) +guides(alpha=FALSE),
       height = 6, width = 8, units = "in")


# Population statistic plots including Manhattan plots -----------------------------------------------------------------
#Plink and vcftools was first run to generate LD and fst values, defining groups according to PCA/DAPC results
#Plink commands to run external to R
#$plink --vcf populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.recode.vcf --double-id --allow-extra-chr --make-bed --out boreogadus.bias
#$plink --bfile boreogadus.bias --allow-extra-chr --r2 --ld-window-kb 100000 --ld-window 1000 --ld-window-r2 0 --out boreogadus.bias

#Remove sites less than 100k bp in proximity, and with a distance greater than 5mbp
boreogadus.LD <- read.csv("boreogadus.bias.LD.100k-5mbp.txt", sep = "\t")
#Average R2 for each SNP using unique SNP identifier
R2.avg <- tapply(boreogadus.LD$R2, boreogadus.LD$SNP_A..CH, mean)
R2.dataset.avg <- as.data.frame(R2.avg)
#Export and link up Chrom and bp information to SNPs
write.csv(R2.dataset.avg, "R2.dataset.avg.100k-5mbp.csv")
R2.dataset.avg <- read.csv("R2.dataset.avg.100k-5mbp.csv")

#Import Fst data
boreogadus.loadings <- read.csv("dapc.dataset.daloadings.outlier.b.grp.csv", sep = ",")
boreogadus.LD <- read.csv("R2.dataset.avg.100k-5mbp.csv", sep = ",")
boreogadus.PCadapt <- read.csv("PCADAPT_k2_results.csv", sep = ",")

#Plot, modifying as needed fot given dataset
gg.boreogadus.pcadapt <- ggplot(boreogadus.loadings, aes(x=Order, y=LD1, col=CHROM)) + 
  geom_point() + theme_bw() + 
  #scale_colour_manual(values=c("#440154", "#471164", "#482173", "#472f7d", "#433e85", "#3e4c8a", "#38588c", "#32648e", "#2d708e", "#297a8e", "#25858e", "#21918c", "#1e9b8a", "#21a685", "#2ab07f", "#3bbb75", "#52c569", "#69cd5b", "#86d549", "#5adb36", "#c2df23", "#e2e418", "#fde725")) + 
  #scale_colour_manual(values=c("#1e9b8aff", "#482173", "#1e9b8aff", "#482173", "#1e9b8aff", "#482173", "#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173", "#1e9b8aff")) + 
  #Colours for b outliers missing chromosomes
  scale_colour_manual(values=c("#482173", "#1e9b8aff", "#482173", "#1e9b8aff", "#482173", "#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173","#1e9b8aff", "#482173", "#1e9b8af5")) + 
  xlab("Chromosome") #+ scale_y_reverse() + geom_hline(yintercept = -5.23, linetype="dashed", col="red") + geom_hline(yintercept = -2.2, linetype="solid", col="red")

gg.boreogadus.pcadapt

#Use facet to plot LD for chromosomes seperately
facet(gg.boreogadus.ld, facet.by = "CHR_A") + guides(colour = FALSE) + ylim(0,0.5)#+ theme(legend.position = c(0.32,0.1))

ggsave(file.path(paste0("boreogadus.ld1.10.pdf")), facet(gg.boreogadus.ld, facet.by = "CHR_A") + guides(colour = FALSE) + ylim(0,0.5),
       width = 12, height = 8, unit = "in")

ggsave(file.path(paste0("boreogadus.loading.b.10.pdf")), gg.boreogadus.pcadapt + guides(colour = FALSE), #+ ylim(0,0.5),
       width = 12, height = 3, unit = "in")

#Plot pairwise fst as a heatmap
#Import data as three column tab delimited file, Pop1, Pop2, Fst
boreogadus.fst <- read.csv("boreogadus_fst_pops_7xii23.txt", sep = "\t")

#Sort by population
boreogadus.fst$Pop1 <- factor(boreogadus.fst$Pop1, levels = c("Western Arctic", "Eastern Arctic", "Hudson Bay Complex", "Newfoundland-Labrador Shelves","Gulf of Saint Lawrence"))
boreogadus.fst$Pop2 <- factor(boreogadus.fst$Pop2, levels = c("Western Arctic", "Eastern Arctic", "Hudson Bay Complex", "Newfoundland-Labrador Shelves","Gulf of Saint Lawrence"))

#Plot Fst using viridis
gg.boreogadus.fst.pops <- ggplot(boreogadus.fst, aes(Pop1, Pop2, fill=Fst.avg)) +
  geom_tile() + scale_fill_viridis_c() + theme_bw()
gg.boreogadus.fst.pops

ggsave(file.path(paste0("boreogadus.fst.pops.pdf")), gg.boreogadus.fst.pops, #+ guides(colour = FALSE), #+ ylim(0,0.5),
       width = 12, height = 8, unit = "in")

#Generated inbreeding coefficient values for populations individually through vcftools
#Plot population statistics as violin plots and geom_points as jittered clouds
boreogadus.het <- read.csv("clusters.het.txt", sep="\t")
boreogadus.het$Ecoregion <- factor(boreogadus.het$Ecoregion, levels = c("Western Arctic", "Eastern Arctic", "Hudson Bay Complex", "Newfoundland-Labrador Shelves","Gulf of Saint Lawrence"))

gg.boreogadus.het <- boreogadus.het %>% #filter(Population != "All") %>% 
  ggplot(aes(x=Cluster, y=F, fill=Cluster)) + 
  geom_violin(alpha=0.3) + theme_bw() + 
  geom_jitter(size=2, position=position_jitterdodge(1.5), aes(col=Cluster)) + 
  scale_fill_manual(values=c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")) +
  scale_colour_manual(values=c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"))
gg.boreogadus.het

ggsave(filename = file.path("boreogadus_het.10x.pdf"), 
       plot = gg.boreogadus.het + ggtitle("boreogadus pop 10x 8311 SNPs"), #+ guides(colour = FALSE),
       height = 3, width = 10, units = "in") 

#Plot length and weight measurements as boxplots
boreogadus.trait <- read.csv("2022_Boreogadus_20230822_6xii23.csv", sep=",")
boreogadus.trait$Ecoregion <- factor(boreogadus.trait$Ecoregion, levels = c("Western Arctic", "Eastern Arctic", "Hudson Bay Complex", "Newfoundland-Labrador Shelves","Gulf of Saint Lawrence"))
boreogadus.trait$Region_echantillonnage <- factor(boreogadus.trait$Region_echantillonnage, levels = c("Beaufort_Sea_ecoregion", "Eastern_Arctic_ecoregion", "0A", "0B","1A", "Hudson_Bay_ecoregion", "WAZ", "2G", "2H", "2J", "3K", "4S", "Saguenay"))

gg.boreogadus.trait <- boreogadus.trait %>%
  ggplot(aes(x=Region_echantillonnage, y=Poids_g, fill=Ecoregion)) + 
  geom_boxplot() + theme_bw() + 
  scale_fill_manual(values=c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")) +
  scale_colour_manual(values=c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"))
gg.boreogadus.trait

ggsave(filename = file.path("boreogadus_weight.10x.pdf"), 
       plot = gg.boreogadus.trait + ggtitle("boreogadus pop 10x 8311 SNPs"), #+ guides(colour = FALSE),
       height = 4, width = 8, units = "in") 

# Admixture plots -----------------------------------------------------------------

#Ran Admixture externally with 5-fold cross validation error.
#Input files created using Plink
#$plink --vcf populations.38131snps.507indwArctogadus.H06.DP.single.final.noArc+hy.0-05.recode.vcf --double-id --allow-extra-chr --make-bed --out boreogadus.pop.10x
#Fix .bim for Admixture
#awk '{$1=0;print $0}' boreogadus.species.10x.bim > boreogadus.species.bim.tmp
#mv boreogadus.species.bim.tmp boreogadus.species.bim

admixture.analysis <- "boreogadus.pop.10x"

# Cross-validation results:

CV.res <- data.frame(k = 1:10,
                     CV = NA,
                     stringsAsFactors = F)

for(i in 1:nrow(CV.res)){
  # Which k
  k <- CV.res[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path(paste0("log",k, ".out")))
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

ggsave(filename = file.path(paste0(admixture.analysis, ".10.CV.png")), 
       plot = gg.CV + ggtitle(paste0(admixture.analysis, " cross validation error")),
       height = 3.5, width = 7, units = "in")   

#Assemble dataset at various levels of k
k <- 6
Q.k2.res <-  read.table(file.path(paste0(admixture.analysis,".",2,".Q")))
Q.k3.res <-  read.table(file.path(paste0(admixture.analysis,".",3,".Q")))
Q.k4.res <-  read.table(file.path(paste0(admixture.analysis,".",4,".Q")))
Q.k5.res <-  read.table(file.path(paste0(admixture.analysis,".",5,".Q")))
Q.k6.res <-  read.table(file.path(paste0(admixture.analysis,".",6,".Q")))

Q.res.boreogadus <- bind_rows(cbind(fam$V1, Q.k6.res, K = 6),
                              cbind(fam$V1, Q.k5.res, K = 5),
                              cbind(fam$V1, Q.k4.res, K = 4),
                              cbind(fam$V1, Q.k3.res, K = 3),
                              cbind(fam$V1, Q.k2.res, K = 2))

head(Q.res.boreogadus)

#Q.fas.res <- cbind(fam.fas$V1, Q.fas.res)

names(Q.res.boreogadus) <- c("ID_GQ", paste0("Q", 1:6), "K")


#Add regional information externally
write.csv(Q.res.boreogadus, paste0(admixture.analysis,".csv"))
Q.res.boreogadus <- read.csv(paste0(admixture.analysis, ".10x.k2-4_7xii23.csv"))

Q.res.boreogadus$Ecoregion <- factor(Q.res.boreogadus$Ecoregion, levels = c("Arctogadus", "Potential hybrid", "Western Arctic", "Eastern Arctic", "Hudson Bay Complex", "Newfoundland-Labrador Shelves","Gulf of Saint Lawrence"))

gg.str.all <- Q.res.boreogadus %>% pivot_longer(cols =  paste0("Q", 1:4), names_to = "Group", values_to = "Q") %>% 
  #  mutate(Group = factor(Group, levels = c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10"))) %>% 
  #  dplyr::filter(ID_GQ %nin% bad.samples) %>% 
  #left_join(pop.data.ddrad) %>% 
  #ggplot(aes(x = reorder(ID_GQ, Longitude_echantillonnage_DD, FUN = function(x) max(as.numeric(x))), y = Q, fill = Group)) + 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  
  geom_col() +
  #facet_grid(. ~Lieu_echantillonnage + Mois_echantillonnage, space = "free", scale = "free")  +
  facet_grid(K ~ Ecoregion, space = "free", scale = "free") +
  scale_fill_manual(values=c("#31688e", "#35b779", "#440154", "#fde725"))
#scale_fill_manual(values=c("#2ed7e8ff", "darkblue"))
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
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") ) + theme(strip.text.x = element_text(size=10))
gg.str.all

ggsave(file.path(paste0("Admixture", admixture.analysis,".10.k2-4.pdf")), plot = gg.str.all + ggtitle(paste0(admixture.analysis, " admixture pop k4")),
       width = 12, height = 6, unit = "in")

# END -----------------------------------------------------------------

