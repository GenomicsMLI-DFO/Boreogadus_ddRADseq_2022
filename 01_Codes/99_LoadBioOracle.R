install.packages("sdmpredictors")


# Load package
library(sdmpredictors)
library(corrplot)

# Explore datasets in the package
list_datasets()

# Explore layers in a dataset
A <- list_layers(datasets = "Bio-ORACLE", terrestrial = FALSE, marine = TRUE, freshwater =  FALSE) 

#Download specific layers to the current directory

# Add a 60 min timeout because download is really slow
options(timeout = max(3600, getOption("timeout")))

env.bioOracle <- load_layers(c("BO22_tempmean_bdmean", #"BO22_tempmin_bdmean", "BO22_tempmax_bdmean", 
                               "BO22_tempmean_ss", #"BO22_tempmin_ss", "BO22_tempmax_ss",  
                               "BO22_salinitymean_bdmean", #"BO22_salinitymin_bdmean", "BO22_salinitymax_bdmean", 
                               "BO22_salinitymean_ss", #"BO22_salinitymin_ss", "BO22_salinitymax_ss", 
                               "BO22_dissoxmean_bdmean",#  "BO22_dissoxmin_bdmean",  "BO22_dissoxmax_bdmean", 
                               "BO22_dissoxmean_ss",#  "BO22_dissoxmin_ss",  "BO22_dissoxmax_ss", 
                               "BO_bathymean"
                               
), datadir = file.path("/media/genyoda/Fast_Storage/Projet/BioOracle/") )


env.bioOracle
names(env.bioOracle)

# Check layer statistics


# Extract -----------------------------------------------------------------

pop.data <- read_csv("00_Data/00_FileInfos/Project_Infos_20230822.csv")
pop.data 

Env.data <- pop.data %>%   
  group_by(Numero_unique_groupe) %>% 
  summarise(Lat = mean(Latitude_echantillonnage_DD),
            Long = mean(Longitude_echantillonnage_DD))


Env.ras <- data.frame(raster::extract(env.bioOracle, Env.data[, c("Long", "Lat")]))

Env.r  <- cor(Env.ras[!is.na(Env.ras[,1]),], method = "pearson")

gg.r.env <- ggcorrplot::ggcorrplot(Env.r,   hc.order = TRUE, type = "lower", lab = TRUE) +
  theme_bw(base_size = 8) + theme(axis.title = element_blank(),
                                  legend.position = "bottom",
                                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~"Spearman's correlations")

gg.r.env 

row.names(Env.ras) <- Env.data$Numero_unique_groupe

pca.env <- vegan::rda(Env.ras[!is.na(Env.ras[,1]),], scale = T)

Env.data <- bind_cols(Env.data, Env.ras)

library(vegan)

res.env <- as.data.frame(scores(pca.env , display="sites", scaling=1)) %>% 
  mutate(Numero_unique_groupe = dimnames(scores(pca.env, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data %>% dplyr::select(Numero_unique_groupe, Region_echantillonnage) %>% distinct())


arrow.env.df <- data.frame(name = scores(pca.env, display="species", choices=1, scaling=1) %>% row.names() ,
                           x0 = 0,
                           y0 = 0,
                           xmax = scores(pca.env, display="species", choices=1, scaling=1) %>% as.vector(),
                           ymax = scores(pca.env, display="species", choices=2, scaling=1) %>% as.vector())  
arrow.factor = 1/4


perc <- round(100*(summary(pca.env)$cont$importance[2, 1:2]), 1)


library(ggrepel)
gg.pca.env <-  res.env %>% 
  ggplot(aes(x =PC1, y = PC2)) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  geom_point(aes(fill = Region_echantillonnage), pch = 21, cex = 2)+
  #scale_fill_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta")
  #)+ 
  
  geom_segment(data = arrow.env.df, aes(x = x0, y = y0, xend = xmax * arrow.factor, yend= ymax * arrow.factor),
               arrow = arrow(), size = 0.5)+
  geom_label_repel(data = arrow.env.df  , aes(label = name, x = xmax *arrow.factor, y = ymax * arrow.factor), size = 2, max.overlaps = 20
  ) +
  #  geom_text_repel(aes(label = Gen_ZONE_FG %>% str_remove("SFA-|NAFO-")), size = 3, max.overlaps = 20
  #  ) +
  labs(x = paste0("PC1 (", perc[1], "%)"), 
       y = paste0("PC2 (", perc[2], "%)")) + 
  #annotate("text", label = paste("R²-adjusted =", RsquareAdj(rdaS.northern)$adj.r.squared %>% as.numeric() %>%  round(digits = 3)),
  #          x = -0.2, y = 0.25, vjust = "inward", hjust = "inward") +
  facet_wrap(~"Principale component analysis") +
  theme_bw(base_size = 8) #+ theme(legend.position = "none")

gg.pca.env

gg.env <- ggpubr::ggarrange(gg.r.env,
                            gg.pca.env,
                            labels = LETTERS,
                            widths = c(2,2))
gg.env


