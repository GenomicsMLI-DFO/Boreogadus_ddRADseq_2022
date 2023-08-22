# Info --------------------------------------------------------------------

# Creation of POP INFO files + barcodes for RAD-seq pipeline with STACKS
#
# Audrey Bourret
# 2023-08-11
#

# Library -----------------------------------------------------------------

library(readxl)
library(tidyverse)
library(BDLG.gestion) #to upload Boreogadus data the first time

# Internal functions
for(i in 1:length( list.files("./01_Codes/Functions") )){
  source(file.path("./01_Codes/Functions",  list.files("./01_Codes/Functions")[i]))
}


# Data --------------------------------------------------------------------

list.files(get.value("info.path"))

IBIS.barcode <-  read_excel(file.path(get.value("info.path"),"BarcodesIBIS.xlsx"))
IBIS.barcode

# Ce qui a été envoyé

GQ.data <- load_DB("Metadata_Analyses_Externes")
GQ.data %>% dplyr::pull(Nom_projet) %>% table()
GQ.data <- GQ.data |> dplyr::filter(Nom_projet == "BOREOGADUS_2022")

#write_csv(GQ.data, "./00_Data/00_FileInfos/2022_Boreogadus_20230811.csv")

# Add infos on the extraction

lab.data <- load_DB("04_Extraits_ADN_ARN")

GQ.data2 <- GQ.data |> left_join(lab.data |> 
                                   dplyr::select("Numero_unique_extrait", "Annee_extraction", "Mois_extraction", "Jour_extraction", "Type_extraction", "Kit_extraction", "Protocole_extraction",    "Responsable_extraction" ,"Notes_extraitADN"))

#write_csv(GQ.data, "./00_Data/00_FileInfos/2022_Boreogadus_20230811.csv")
write_csv(GQ.data2, "./00_Data/00_FileInfos/2022_Boreogadus_20230822.csv")

GQ.data <- readr::read_csv("./00_Data/00_FileInfos/2022_Boreogadus_20230822.csv")

GQ.data |> group_by(No_soumission_GQ) |> group_by(No_plaque_envoi) |> summarise(N = n())

# Identify the duplicated samples
names(GQ.data)
table(duplicated(GQ.data$Numero_unique_specimen))
table(duplicated(GQ.data$Numero_unique_extrait))


GQ.data <- GQ.data %>% mutate(Dup_specimen = duplicated(Numero_unique_specimen),
                              Cat_sample = ifelse( Dup_specimen == T , "Duplicate", "Sample"),
                              ID_GQ = ifelse(Cat_sample == "Duplicate", paste0(Numero_unique_specimen, "_rep"), Numero_unique_specimen),
                              ID_GQ = ifelse(duplicated(ID_GQ), paste0(ID_GQ, "_1"), ID_GQ)) |>
  left_join(IBIS.barcode %>% dplyr::select(No_puits_envoi = Cell2, Barcode))



table(duplicated(GQ.data$ID_GQ))
table(GQ.data$Cat_sample)

GQ.data$ID_GQ[duplicated(GQ.data$ID_GQ)]

#GQ.data %>% View()

GQ.data %>%  dplyr::mutate(puit_range = stringr::str_sub(No_puits_envoi, 1,1),
                           puit_col = stringr::str_sub(No_puits_envoi, 2,3)) |>
  ggplot(aes(x = puit_col, y = puit_range, fill = Cat_sample)) +
  geom_bin2d() +
  ggtitle("Duplicated on plate") +
  scale_y_discrete(limits=rev) +
  facet_wrap(~No_plaque_envoi) +
  theme_bw()

write_csv(GQ.data, file = file.path(get.value("info.path"),"Project_Infos_20230822.csv"))

# Create the first popmap

pop.info <- GQ.data %>% select(ID_GQ) %>% mutate(POP = "NoPop")
pop.info

write.table(pop.info,
            file = file.path(get.value("info.path"), "popmap.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)


# Create barcode files
## Check that plate ID are the same as the files
GQ.data %>% pull(No_plaque_envoi) %>% unique()
## Check that every individual as an unique ID

GQ.data %>% pull(ID_GQ) %>% length()
GQ.data %>% pull(ID_GQ) %>% unique() %>% length()
GQ.data %>% filter(duplicated(ID_GQ))  # pull(ID_GQ[duplicated(ID_GQ)])# %>% duplicated()

# Loop to create all the files relating individuals to barcodes
for(x in unique(GQ.data$No_plaque_envoi)){

   write.table(GQ.data %>% filter(No_plaque_envoi == x) %>%
                              select(Barcode, ID_GQ),
                              file = file.path(get.value("info.path"), paste(x,"barcodes.txt", sep = "_")),
                              quote = FALSE, sep = "\t",
                             row.names = F, col.names = F)

}

