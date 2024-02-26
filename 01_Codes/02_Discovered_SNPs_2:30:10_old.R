# Info --------------------------------------------------------------------

# RAD-seq pipeline with STACKS version >2, up-to Gstacks
# NS data from 1 batch
# Trimmomatics + no rad check on one side (cut 5 pb) MAYBE
# With and without ref genome
#
# Audrey Bourret
# 2023-08-14
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

# FastQC ------------------------------------------------------------------

# Took around 30 min by plate (NS and HI)
# Now with an overwrite function

fastqc <- function(folder.in, folder.out, overwrite = F, nthread) {
  #@folder.in : where to find the data
  #@folder.out : where to put the data
  #@overwrite : if T, will remove everything within the folder  
  
  if(get_os() %in% c("os","linux")){ # to run only on os and not windows
    
    cat("Performing a FastQC analysis\n")
    
    # Remove old files or not (addition June 2020)
    if(isTRUE(overwrite)){
      file.remove(list.files(folder.out, full.name = T, pattern ="fastqc"))
      files.to.use <- list.files(folder.in, full.names = T, pattern = "fastq") %>% 
        str_remove(".md5") %>% 
        str_subset("R1.fastq|R2.fastq") %>% unique()
    }
    
    if(isTRUE(overwrite == F)){
      previous.analysis <- list.files(folder.out, pattern = ".html") %>% str_replace("_fastqc.html", ".fastq.gz")
      
      files.to.use <- list.files(folder.in, full.names = T, pattern = "fastq") %>% 
        str_remove(".md5") %>% unique() %>% 
        str_subset(paste(previous.analysis, collapse = "|"), negate = T) %>% 
        str_subset("R1.fastq|R2.fastq")
      
    }
    
    cat("Results could be find here:", folder.out ,"\n")
    
    mclapply(files.to.use,
             FUN = function(x){
               
               cmd <- paste("--outdir", folder.out, x, 
                            "-t", 1)
               system2("fastqc", cmd) 
               
             } ,
             mc.cores = nthread
    )
    
    
  } else {cat("Cannot perform FastQC on windows yet -- sorry!!")}
} # End of my function

# Test a multiqc function

multiqc <- function(folder.in){
  # Multi QC aggregation - run pretty fast ...
  for(s in c("R1", "R2")){
    print(s)  
    cmd <- paste(paste(list.files(folder.in, full.names = T) %>% 
                         #str_subset(l) %>% 
                         str_subset(paste0("_",s)) %>% 
                         str_subset(".zip"), collapse = " "),
                 "--outdir", file.path(folder.in, "MultiQC_report"),
                 "--filename", paste0("multiqc_report_", s, ".html"),
                 "-f" # to rewrite on previous data
    )
    
    system2("multiqc", cmd)
    
  } 
  
}



# Run FastQC - change for overwrite T for the first time to or it will not work
fastqc(folder.in = get.value("raw.path"), folder.out = get.value("fastqc.raw.path"), overwrite = F, nthread = 20)


# Multi QC aggregation - run pretty fast ...

multiqc(folder.in = get.value("fastqc.raw.path"))

# CHECK on the R2 side that the adapter "CGG" is present. 
# IF not, the PART2 of trimmomatic allows to remove 5 pb on this side

# Trimmomatics ------------------------------------------------------------

# In paired-end
# To remove the Illumina adapter and to cut 3pb in R2

# Check that you can reach the program
trimmomatic.path <- "/home/genyoda/Documents/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar"
system2("java", paste("-jar", trimmomatic.path, "PE", "-version"), stdout=T, stderr=T) 

# Check which files you want to use

files.to.use <- list.files(get.value("raw.path"), full.names = T) %>% 
  str_subset("R1") %>% 
  str_subset(".md5", negate = T)

#start_time3 <- Sys.time()

### BENCHMARK ###
# On 4 files of 1000 reads
# 80 cores in Trimmomatics: 1.249837 secs
# 8 cores in Trimmomatics: 1.179123 secs 
# 4 cores in Trimmomatics / parallele in R 4 cores: 0.3734035 secs
#  1 core in Trimmomatics / parallele in R 4 cores: 0.3527803 secs

# With 60 cores, memory is full, so 20 seems a better compromise
# checked in Terminal with "free -m" or "top"

# # STEP 1 - REMOVE 3 PB on R2 (must do it first because of misalignment)
# 
# mclapply(files.to.use,
#          FUN = function(x){
#            
#            # Names of the important files
#            file1 <- x
#            file2 <- file1 %>% str_replace("_R1", "_R2")
#            
#            file2.out <- file2 %>% str_replace(get.value("raw.path"), get.value("trimmo.path")) %>% 
#              str_replace("_R2.fastq.gz", "_R2_HC3.fastq.gz")
#            
#            logout <- file1 %>% str_replace(get.value("raw.path"), get.value("trimmo.log")) %>% 
#              str_replace("_R1.fastq.gz", ".log")
#            
#            # The command
#            cmd1 <- paste("-jar",
#                          trimmomatic.path, "SE",
#                          "-threads", 1,
#                          #"-trimlog", logout %>% str_replace(".log", "_HC.log"),
#                          file2, file2.out,
#                          #"ILLUMINACLIP:00_Data/00_FileInfos/adapters/TruSeq3-PE-2.fa:2:30:10",
#                          "HEADCROP:3",
#                          sep = " ")
#            
#            A1 <- system2("java", cmd1, stdout=T, stderr=T) # to R console
#            A1
#            
#            # save a file log
#            cat(file = logout,
#                #"STEP1 - REMOVE ILLUMINA ADAPTORS", cmd1, A1,
#                "STEP1 - CUT 3 pb ON R2 PAIRED because of quality drop on Novaseq",cmd1, A1, # what to put in my file
#                append= F, sep = "\n\n")
#            
#          } ,
#          mc.cores = 20
# )
# 
# 

# STEP 2 - REMOVE Illumina adaptors

mclapply(files.to.use,
         FUN = function(x){
           
           # Names of the important files
           file1 <- x
           file2 <- file1 %>% str_replace("_R1", "_R2")
           
           fileout <- file1 %>%  str_replace(get.value("raw.path"), get.value("trimmo.path")) %>% 
             str_replace("_R1.fastq.gz", ".fastq.gz")   
           
           logout <- file1 %>% str_replace(get.value("raw.path"), get.value("trimmo.log")) %>% 
             str_replace("_R1.fastq.gz", ".log")
           
           # The command
           
           cmd1 <- paste("-jar",
                         trimmomatic.path, "PE",
                         "-threads", 1,
                         #"-trimlog", logout,
                         file1, file2,
                         "-baseout", fileout,
                         "ILLUMINACLIP:00_Data/00_FileInfos/adapters/TruSeq3-PE-2.fa:2:30:10",
                         #"HEADCROP:1-99",
                         sep = " ")
           
           A1 <- system2("java", cmd1, stdout=T, stderr=T) # to R console
           A1
           
           # save a file log
           cat(file = logout,
               "STEP2 - REMOVE ILLUMINA ADAPTORS", cmd1, A1,
               #"STEP1 - CUT 3 pb ON R2 PAIRED because of quality drop on Novaseq",cmd1, A1, # what to put in my file
               append= T, sep = "\n\n")
           
         } ,
         mc.cores = 20
)


# SOME POST ANALYSIS FILE MANIPULATION

# Rename the files

old.name <-list.files(get.value("trimmo.path"), full.names = T)
old.name

new.name <- old.name %>% 
  str_replace("_1P.fastq", "_R1.fastq") %>% 
  str_replace("_2P.fastq", "_R2.fastq") 

new.name

file.rename(from = old.name,
            to = new.name)

# Then remove ALL unecessary files

file.to.remove <- list.files(get.value("trimmo.path"), full.names = T, pattern = "1U.fastq.gz|2U.fastq.gz|HC3.fastq.gz")
sum(file.size(file.to.remove))/  1024 / 1024 / 1024
file.to.remove

file.remove(file.to.remove)


# Compute stats by plate/RunSeq


# Will work for paired reads

trim.data <- data.frame(ID = character(),
                        RunSeq = character(),
                        No_plaque = character(),
                        Ntotal = numeric(),
                        Nsurvival = numeric(),
                        stringsAsFactors = F)

for(x in  list.files(get.value("trimmo.log"), pattern = ".log")){
  
  ID <- x %>% str_remove(".log")
  ID.int = ID %>% str_replace("[.][A-Z][:digit:][:digit:][:digit:][.]|[.][A-Z][:digit:][:digit:][:digit:]---[A-Z][:digit:][:digit:][:digit:][.]", "___")
  RunSeq = sapply(str_split(ID.int, "___"), `[`, 1)
  No_plaque = sapply(str_split(ID.int, "___"), `[`, 2)
  temp <- readLines(file.path(get.value("trimmo.log"), x ))
  temp <- temp %>% str_subset(pattern = "Input Read Pairs")
  Ntotal    <- sapply(str_split(temp, " "), `[`, 4)
  Nsurvival <- sapply(str_split(temp, " "), `[`, 7)
  
  trim.data <- bind_rows(trim.data,
                         data.frame(ID = ID,
                                    RunSeq = RunSeq,
                                    No_plaque = No_plaque,
                                    Ntotal = as.numeric(Ntotal),
                                    Nsurvival = as.numeric(Nsurvival),
                                    stringsAsFactors = F))
  
}
head(trim.data)
trim.data <- trim.data %>% filter(Ntotal != Nsurvival)

write_csv(x = trim.data, file = file.path(get.value("trimmo.log"), "Summary_Nreads.csv") )

gg.trim <- trim.data %>% ggplot(aes(x = No_plaque, y = Nsurvival, fill = RunSeq)) +
  geom_bar(stat = "identity") +
  labs(title = "N reads by plates post trimmomatic") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

gg.trim <- trim.data %>% mutate(perc_surv = Nsurvival/Ntotal) %>% 
  ggplot(aes(x = perc_surv, fill = RunSeq)) +
  geom_histogram() +
  labs(title = "N reads by plates post trimmomatic") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

gg.trim <- trim.data %>% mutate(perc_surv = Nsurvival/Ntotal) %>% 
  ggplot(aes(x = Nsurvival, fill = RunSeq)) +
  geom_histogram() +
  labs(title = "N reads by plates post trimmomatic") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

gg.trim <- trim.data %>% mutate(perc_surv = Nsurvival/Ntotal) %>% 
  ggplot(aes(x = Ntotal, y = Nsurvival, col = RunSeq)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(title = "N reads by plates post trimmomatic") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

ggsave(filename = file.path(get.value("trimmo.log"), "Summary_Nreads.png"),
       plot = gg.trim,
       width = 5, height = 4, units = "in")


# If you want the log file to be ignore, run the following :

cat("*.log", "!.gitignore", sep = "\n",
    file = file.path(get.value("trimmo.log"), ".gitignore"))


# Run FastQC - change for overwrite T first the first time to or it will not work
# THIS PART SHOULD BE BENCHMARKED TOO

fastqc(folder.in = get.value("trimmo.path"), folder.out = get.value("fastqc.trim.path"), 
       overwrite = T, nthread = 20)

multiqc(get.value("fastqc.trim.path"))


# Demultiplex -------------------------------------------------------------

# Barcode files are important for thise step
# They couldbe created with the 01_Create_PopMap.R

# TAKE CARE, all individuals with similar names will be collapse
length(pop.info$Sample)
length(pop.info$Sample %>% unique())

# If you specify the same output directory for two differents
# process_radtags runs, the second run will overwrite identical filenames
# from the first run. Output the data into separate directories, and then
# concatenate the shared samples together after the runs complete

sub.dir <- c("NS.2050.002")
#sub.dir <- "NS.1795.003"
sub.dir[3]

# Before running, check that the right barcode file is found - here an example
# not in a loop

files.to.use <-  list.files(get.value("trimmo.path"), full.names = T, pattern = "_R1") #%>% 
#  str_subset(sub.dir[1])
files.to.use

barcode <- list.files(get.value("info.path"), full.names = T, pattern = "barcodes.txt") %>% 
  str_subset(files.to.use[1] %>% str_remove(get.value("trimmo.path")) %>% 
               str_remove("_R1.fastq.gz") %>% 
               str_remove(sub.dir[1]) %>%
              #str_remove("[.][A-Z][:digit:][:digit:][:digit:][.]") %>% 
               str_remove("[.][A-Z][:digit:][:digit:][:digit:]---[A-Z][:digit:][:digit:][:digit:][.]") %>% 
               str_remove("/") %>% 
               str_replace("Pradseq_2022", "Pradseq_22")) 
barcode



for(i in sub.dir){
  
  files.to.use <-  list.files(get.value("trimmo.path"), full.names = T, pattern = "_R1") %>% 
    str_subset(i) #%>% str_subset("P02|P03")
  
  # files.to.use
  
  cat(paste("\nWorking with the run:", i),
      paste("There are" , length(files.to.use) , "plates within this run"),
      sep= "\n")
  
  # Créer un sous-dossier s'il n'existe pas
  if(file.exists(file.path(get.value("demulti.path"), i)) == F) {
    
    cat(paste("\nCreating a new directory:", file.path(get.value("demulti.path"), i)),
        sep= "\n")
    
    dir.create(file.path(get.value("demulti.path"), i), recursive = T)}
  
  
  # Parallel version of process_radtag 
  
  
  mclapply(files.to.use,
           FUN = function(x){
             
             # Files
             file1 <- x
             file2 <- file1 %>% str_replace("_R1", "_R2")
             
             barcode <- list.files(get.value("info.path"), full.names = T, pattern = "barcodes.txt") %>% 
               str_subset(x %>% str_remove(get.value("trimmo.path")) %>% 
                            str_remove("_R1.fastq.gz") %>% 
                            str_remove(i) %>%
                            #str_remove("[.][A-Z][:digit:][:digit:][:digit:][.]") %>% 
                            str_remove("[.][A-Z][:digit:][:digit:][:digit:]---[A-Z][:digit:][:digit:][:digit:][.]") %>% 
                            str_remove("/")%>% 
                            str_replace("Pradseq_2022", "Pradseq_22")) 
             
             plate <- barcode %>% str_remove(get.value("info.path")) %>% 
               str_remove("_barcodes.txt") %>% 
               str_remove("/")
             
             # Créer un sous-dossier s'il n'existe pas
             if(file.exists(file.path(get.value("demulti.path"), i, plate)) == F) {
               
               #cat(paste("\nCreating a new directory:", file.path(get.value("demulti.path"), i)),
               #     sep= "\n")
               
               dir.create(file.path(get.value("demulti.path"), i, plate), recursive = T)}
             
             # Command
             cmd <- paste("--paired",
                          "-1", file1,
                          "-2", file2,
                          "-o", file.path(get.value("demulti.path"), i, plate),
                          "--inline_null",   
                          "-b", barcode,
                          "--renz_1", "pstI", # Check RAD site on R1 
                          "--renz_2", "mspI", # CGG on R2 
                          "-E", "phred33",
                          "--filter-illumina",
                          "-c", # clean
                          "-r", #rescue
                          "-q", #check quality
                          "-t", 135,  #truncate at 150 - 6 (restriction site) - 8 (max barcode) = 136 (all reads within sample must be the same length)
                          "-i", "gzfastq"
             )
             
             A <- system2("process_radtags", cmd, stdout=T, stderr=T)
             A
             # save a log file 
             log.file <- file1 %>% str_replace(get.value("trimmo.path"), get.value("demulti.log")) %>% 
               str_replace(".fastq.gz","_summary.log") %>% 
               str_remove("_R1")
             
             cat(file = log.file,
                 cmd, "\n\n",
                 A, # what to put in my file
                 append= F, sep = "\n")
             
             # Detailed log file
             file.rename(from = list.files(file.path(get.value("demulti.path"), i, plate), full.names = T, pattern = ".log"),
                         to = log.file %>% str_replace("summary", "detailed")
             )
             
           } ,
           mc.cores = 10
  )
  
  gc()
  
}


# If you want the log/csv file to be ignore, run the following :

cat("*", "!.gitignore", "!*.png", "!AllIndividuals_Nreads.csv", sep = "\n",
    file = file.path(get.value("demulti.log"), ".gitignore"))



# Extract data for each summary_detailed
# You must change the str_subset code that is the index for each project

for(x in list.files( get.value("demulti.log"), pattern = "detailed", full.names = T)){
  
  data <- readLines(x) %>% 
    str_subset("S_") %>% 
    str_split(pattern = "\t")
  
  cat("\n",length(data), "samples retrieved in", x %>% str_remove(get.value("demulti.log")) %>% str_remove("/"))
  
  data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
  # If there is a RUN column, it's because some individuals come from more than one sample
  names(data) <- c("Barcode", "Filename", "Total", "NoRadTag", "LowQuality", "Retained", "Pct_Retained", "Pct_Total_Reads")
  
  data <- data %>% mutate(Total = as.numeric(as.character(Total)),
                          NoRadTag = as.numeric(as.character(NoRadTag)),
                          LowQuality = as.numeric(as.character(LowQuality)),
                          Retained = as.numeric(as.character(Retained)),
                          Run = x %>% str_remove(get.value("demulti.log")) %>% 
                            str_remove("_detailed.log")  %>% 
                            str_remove("/")
  )
  write_csv(data, x %>% str_replace("detailed.log", "Nreads.csv"))
  
}

# Create one big log file 

Nreads.data <- data.frame()

for(x in list.files( get.value("demulti.log"), pattern = "Nreads", full.names = T) %>% str_subset(".csv") %>% str_subset("NS.")){
  data.int <- read_csv(x)
  
  Nreads.data <- bind_rows(Nreads.data, data.int)
  
}

# Compute N read removed
Nreads.data$Removed <- Nreads.data$Total - Nreads.data$Retained  

# Add a column for the run name
Nreads.data$RunSeq <- paste(sapply(str_split(Nreads.data$Run, "[.]"),`[`,1),
                            sapply(str_split(Nreads.data$Run, "[.]"),`[`,2),
                            sapply(str_split(Nreads.data$Run, "[.]"),`[`,3),
                            sep = ".")


# Save the result
#write_csv(Nreads.data, file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))

Nreads.data <- read_csv(file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))

trim.data <- read_csv(file = file.path(get.value("trimmo.log"), "Summary_Nreads.csv") ) %>% 
           mutate( No_plaque =  No_plaque %>% str_replace("Pradseq_2022", "Pradseq_22"))

Nreads.data <- Nreads.data %>% left_join(pop.data, by = c("Filename" = "ID_GQ")) 

# Check that we have all the data - all with unique data

head(Nreads.data)
Nreads.data %>% nrow() / length(sub.dir)
Nreads.data$Filename %>% unique() %>% length()

# Stats by plate
Nreads.data %>% group_by(No_plaque_envoi) %>% summarise(Retained = sum(Retained)) %>%
  group_by() %>% 
  summarise(mean = mean(Retained)/2,
            sd = sd(Retained)/2,
            min = min(Retained)/2,
            max = max(Retained)/2)

# Impact of overall process

gg.radtag <- Nreads.data %>% group_by(No_plaque_envoi, RunSeq) %>% 
  #mutate(RunSeq = paste()) 
  summarise(Total = sum(Total)/2,
            Retained = sum(Retained)/2) %>% 
  left_join(trim.data  %>%  select(RunSeq, No_plaque_envoi =  No_plaque, Nsurvival)) %>% 
  mutate(perc_with_barcode = Total/Nsurvival,
         perc_keep = Retained / Nsurvival,
         perc_keep_over_barcode = Retained / Total) %>% 
  select(RunSeq, No_plaque_envoi, perc_with_barcode, perc_keep,  perc_keep_over_barcode) %>% 
  pivot_longer(c(perc_with_barcode, perc_keep,  perc_keep_over_barcode), names_to = "Cat", values_to = "Perc") %>% 
  ggplot(aes(x = No_plaque_envoi, y = Perc, col = RunSeq, group = RunSeq)) +
  geom_point(position= position_dodge(width = 1/4)) +
  labs(title = "% reads retained post demultiplexing") +
  facet_grid(Cat ~ . ) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.radtag

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Preads_byPlates.png"),
       plot = gg.radtag,
       width = 8, height = 6, units = "in")


gg.retained.plate <- Nreads.data %>% #gather(Removed, Retained, key = "Cat", value = "N") %>% #head()
  group_by(RunSeq, No_plaque_envoi) %>% summarise(Retained = sum(Retained)/2,
                                                  Removed = sum(Removed)/2) %>% 
  ggplot(aes(y = Retained, x = No_plaque_envoi, fill = RunSeq)) + 
  geom_bar(stat = "identity")+
  theme_bw() + 
  labs(title = "N reads retained post demultiplexing") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")
gg.retained.plate

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byPlates.png"),
       plot = gg.retained.plate,
       width = 5, height = 4, units = "in")


# Check that we have all the data - all with unique data

head(Nreads.data)
Nreads.data %>% nrow() / length(sub.dir)
Nreads.data$Filename %>% unique() %>% length()

# Stats by IND
Nreads.data %>% group_by(Filename) %>% summarise(Retained = sum(Retained)) %>%
  group_by() %>% 
  summarise(mean = mean(Retained)/2,
            sd = sd(Retained)/2,
            min = min(Retained)/2,
            max = max(Retained)/2)


gg.retained.ind <- Nreads.data %>% gather(Removed, Retained, key = "Cat", value = "N") %>% #head()
  group_by(Filename, No_plaque_envoi, Cat) %>% summarise(N = sum(N)/2) %>% 
  ggplot(aes(x = Filename, y = N, fill= Cat)) + 
  geom_bar(stat= "identity") +
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000),
  #                   labels = c(1:7))+
  labs(y = "N reads", x = "Samples", title = "N reads by ind") +
  facet_wrap(~ No_plaque_envoi, scale = "free_x") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        legend.position = "bottom")
gg.retained.ind

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byInd.png"),
       plot = gg.retained.ind,
       width = 8, height = 8, units = "in")


gg.retained.ind2.a <- Nreads.data %>% 
  mutate(barcode_length = nchar(Barcode.y),
         first.nuc = str_sub(Barcode.y, 1,1)) %>% 
  group_by(Filename, barcode_length, first.nuc) %>% summarise(Retained = sum(Retained)/2) %>%
  #  arrange(Retained) %>% 
  ggplot(aes(x = reorder(Filename, Retained), y = Retained)) +
  scale_y_continuous(trans = "log10") +
  geom_bar(stat = "identity") +
  #facet_grid(barcode_length ~ first.nuc, scale = "free_x", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank())
gg.retained.ind2.a



gg.retained.ind2.b <- Nreads.data %>% 
  mutate(barcode_length = nchar(Barcode.y),
         first.nuc = str_sub(Barcode.y, 1,1)) %>% 
  group_by(Filename, barcode_length, first.nuc) %>% summarise(Retained = sum(Retained)/2) %>%
  #  arrange(Retained) %>% 
  ggplot(aes(x = reorder(Filename, Retained), y = Retained )) +
  scale_y_continuous(trans = "log10") +
  geom_bar(stat = "identity") +
  facet_grid(barcode_length ~ first.nuc, scale = "free_x", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank())
gg.retained.ind2.b


gg.retained.ind2 <- ggpubr::ggarrange(gg.retained.ind2.a, gg.retained.ind2.b,
                                      nrow = 2, heights = c(1:3))
gg.retained.ind2

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byInd_withBarcodes.png"),
       plot = gg.retained.ind2,
       width = 8, height = 8, units = "in")

# gg.retained.pop <- Nreads.data %>% 
#   group_by(Filename, Objective, Pop, No_plaque_envoi) %>% summarise(Retained = sum(Retained)/2) %>% 
#    ggplot(aes(x = Pop, y = Retained)) + 
#  geom_hline(yintercept = 0) +
#     geom_boxplot() +
#   geom_jitter(aes(col= factor(No_plaque_envoi)), cex = 1, alpha = 0.5) +
#   labs(y = "N reads (log)") +
#   scale_y_continuous(trans = "log10") +
#   facet_grid(.~Objective, space = "free", scales = "free") +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         legend.position = "none")
# 
# gg.retained.pop
# 
# #ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byPop.png"),
# #       plot = gg.retained.pop,
# #       width =7, height = 5, units = "in")

gg.barcodes <- Nreads.data %>% 
  group_by(No_plaque_envoi, Barcode.y, Filename,RunSeq) %>% summarise(Total = sum(Total)/2) %>% 
  ggplot(aes(x = Barcode.y, y = Total)) + 
  geom_boxplot() +
  geom_jitter(aes(col =RunSeq), cex = 1, alpha = 0.5) +
  # facet_wrap(~Barcode.y)
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000),
  #                   labels = c(1:10))+
  #  scale_y_continuous(trans = "log10") +
  labs(y = "N total reads (log)", x = "Barcodes") +
  # facet_grid(. ~ Espece, scale = "free_x", space = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 ),
        legend.position = "bottom")
gg.barcodes

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byBarcode.png"),
       plot = gg.barcodes,
       width =8, height = 5, units = "in")

# This is THE graph
gg.barcodes.summary <- Nreads.data %>% 
  group_by(Barcode.y, Filename) %>% summarise(Total =mean(Total)/2, 
                                              Retained =mean(Retained)/2,
                                              NoRadTag = mean(NoRadTag)/2,
                                              LowQuality = mean(LowQuality)/2) %>%
  pivot_longer(c(Total, Retained, NoRadTag, LowQuality), names_to = "Cat", values_to = "N") %>% 
  mutate(barcode_length = nchar(Barcode.y),
         first.nuc = str_sub(Barcode.y, 1,1)) %>% 
  ggplot(aes(x = barcode_length, y =N, group = barcode_length) )+ 
  geom_violin() +
  geom_jitter(alpha = 0.5, aes(col = first.nuc), cex = 1) +
  #geom_violin()+
  # facet_wrap(~Barcode.y)
  scale_y_continuous(trans = "log10") +
  labs(y = "N reads (log)", x = "Barcodes length (4-8)") +
  facet_grid(Cat ~ ., scale = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 ),
        legend.position = "bottom")

gg.barcodes.summary

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byBarcode_Summary.png"),
       plot = gg.barcodes.summary,
       width =6, height = 6, units = "in")

gg.barcodes.summary2 <- Nreads.data %>% 
  group_by(RunSeq,Barcode.y, Filename) %>% summarise(Total =mean(Total)/2, 
                                                     Retained =mean(Retained)/2,
                                                     NoRadTag = mean(NoRadTag)/2,
                                                     LowQuality = mean(LowQuality)/2) %>%
  pivot_longer(c(Total, Retained, NoRadTag, LowQuality), names_to = "Cat", values_to = "N") %>% 
  mutate(barcode_length = nchar(Barcode.y),
         first.nuc = str_sub(Barcode.y, 1,1)) %>% 
  group_by(RunSeq, first.nuc, barcode_length) %>% summarise(Nretained = mean(N[Cat == "Retained"])) %>% 
  arrange(Nretained) %>% 
  ggplot(aes(x = factor(barcode_length), y = first.nuc, fill = Nretained)) +
  geom_bin2d() + theme_bw() + facet_wrap(~RunSeq)
gg.barcodes.summary2

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byBarcode_Summary2.png"),
       plot = gg.barcodes.summary2,
       width =5, height = 4, units = "in")#


## Create one file with up to 3 files by individual
## UPDATE: sould work with up to 3 subdir 
sub.dir <- c("NS.2050.002")


sub.dir
plates <- file.path(get.value("demulti.path"), sub.dir[1]) %>% list.files()
plates

for(p in plates){
  
  files.to.use <- list.files(file.path(get.value("demulti.path"), sub.dir[1], p), full.names = T)
  
  cat(paste("\nWorking with the plate:", p),
      paste("There are" , length(files.to.use) , "files to process"),
      sep= "\n")
  
  mclapply(files.to.use,
           FUN = function(x){
             
             # Files
             file1 <-  x
             
             if(length(sub.dir) >= 2){
               file2 <- file1 %>% str_replace(sub.dir[1], sub.dir[2])  
             }
             
             if(length(sub.dir) >= 3){
               file3 <- file1 %>% str_replace(sub.dir[1], sub.dir[3])  
             }
             
             file.join <- file1 %>% str_remove(file.path(sub.dir[1], p))# %>% 
             #str_remove("[:digit:][:digit:][:digit:][:digit:][.][:digit:][:digit:][:digit:][.][A-Z][:digit:][:digit:][:digit:][.]Panomics-" )
             
             if(length(sub.dir) == 1){
               
              file.rename(from = file1, to = file.join) 
               
             }
             
             if(length(sub.dir) == 2){
               # Cat command - this is so simple !
               cmd <- paste(file1, file2, ">", file.join)
               system2("cat", cmd, stdout=T, stderr=T)
             }
             
             if(length(sub.dir) == 3){
               # Cat command - this is so simple !
               cmd <- paste(file1, file2, file3, ">", file.join)
               system2("cat", cmd, stdout=T, stderr=T)
             }
             
             
           } ,
           mc.cores = 20
  )
  
  gc()
  
}

# Number of files observed
list.files(file.path(get.value("demulti.path")), pattern = ".fq.gz") %>% length()

# Number of files expected
384 * length(plates)

# FastQC

#fastqc(folder.in = get.value("demulti.path"), folder.out = get.value("fastqc.demulti.path"), 
#       overwrite = T, nthread = 20)

#multiqc(get.value("fastqc.demulti.path"))


# Select individuals for tests --------------------------------------------

# Select between 10-20 representative individuals
# Will be used both for testing alignement and stacks

Nreads.data <- read_csv(file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))
Nreads.data <- Nreads.data %>% left_join(pop.data, by = c("Filename" = "ID_GQ")) 

# Create the list of individuals for testing parameters :

hist(Nreads.data$Retained)
#summary(Nreads.data$Retained.by.sample)
summary(Nreads.data)
Nreads.data %>% pull(Cat_sample) %>% unique()


Nreads.data %>% group_by(Region_echantillonnage) %>% 
  summarise(N = n()) %>% View()

test.ID <- Nreads.data %>% mutate(RunID = sapply(str_split(Run, "[.]"), `[`, 2)) %>% 
  filter(
    Cat_sample == "Sample"
  ) %>%  
  group_by(Filename, Region_echantillonnage) %>% 
  summarise(Retained = sum(Retained)) %>% 
  filter(Retained >= quantile(Retained, probs = 0.25),
         Retained <= quantile(Retained, probs = 0.75)) %>% 
  group_by(Region_echantillonnage) %>% 
  sample_n(2) %>% pull(Filename)  

test.ID 

write.table(pop.info %>% select(Sample) %>% 
              filter(Sample %in% test.ID) %>% 
              mutate(POP = "noPOP"), 
            file = file.path(get.value("info.path"), "popmap.test_samples.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)


# BWA - index the reference genome ----------------------------------------

# Try first with an old Boreogadus assembly
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_900302515.1/

# Check that 
A <- system2("bwa", "", stdout=T, stderr=T)

list.files(get.value("ref.genome.path"))

cmd <- paste("index",
             "-p",  file.path(get.value("ref.genome.path"),"Bsaida.Genome.August2023"), 
             "-a", "bwtsw",
             file.path("00_Data/99_REF_Genome/ASM90030251v1/GCA_900302515.1/GCA_900302515.1_ASM90030251v1_genomic.fna")
)


A <- system2("bwa", cmd, stdout=T, stderr=T)
A
# save a log file 

cat(file = file.path(get.value("ref.genome.path"), "Bsaida.Genome.August2023.log" ),
    cmd, "\n\n",
    A, # what to put in my file
    append= T, sep = "\n")



# Then, second, garMor3
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902167405.1/

# Check that 
A <- system2("bwa", "", stdout=T, stderr=T)

list.files(get.value("ref.genome.path"))

cmd <- paste("index",
             "-p",  file.path(get.value("ref.genome.path"),"gadMor3.Genome.August2023"), 
             "-a", "bwtsw",
             file.path("00_Data/99_REF_Genome/gadMor3.0/GCA_902167405.1/GCA_902167405.1_gadMor3.0_genomic.fna")
)


A <- system2("bwa", cmd, stdout=T, stderr=T)
A
# save a log file 

cat(file = file.path(get.value("ref.genome.path"), "gadMor3.August2023.log" ),
    cmd, "\n\n",
    A, # what to put in my file
    append= T, sep = "\n")



# BWA - Align reads to the reference genome -------------------------------

# This part take times, can start on a subset (and test Stack at the same times)
# 2-3 min by file

TEST.ID <- read.table( file.path(get.value("info.path"), "popmap.test_samples.txt")) %>% pull(V1) %>% paste(collapse = "|")
#TEST.ID <-  test.ID %>% paste(collapse = "|")

# ATTENTION - ICI POUR PAIRED 
# Test version
demulti.files <- list.files(get.value("demulti.path"), pattern = ".1.fq.gz", full.names = T) %>%
  str_subset(".rem.", negate = T) %>%  
  str_subset(TEST.ID)
# Complete version
demulti.files <- list.files(get.value("demulti.path"), pattern = ".1.fq.gz", full.names = T) %>% 
  str_subset(".rem.", negate = T) #%>%  str_subset(TEST.ID)

demulti.files %>% length()

#demulti.files[1:10]

mclapply(demulti.files,
         FUN = function(x){
           # How I will rename all this : 
           file.R1 <- x
           file.R2 <- x %>% str_replace(".1.fq.gz", ".2.fq.gz")
           file.bam <- x %>% str_replace(get.value("demulti.path"), get.value("align.path")) %>% 
             str_replace(".1.fq.gz", ".bam")
           file.sort.bam <- file.bam %>% str_replace(".bam", ".sorted.bam")
           stat.tsv <- file.bam %>% str_replace(".bam", ".stat.tsv")
           
           # DO THE ALIGMENT  
           if(!file.exists(stat.tsv)){ # do it only if it doesn't exist
             
             cmd1 <- paste("mem",
                           "-t", 16,
                           "-M",
                           
                           file.path(get.value("ref.genome.path"),"gadMor3.Genome.August2023"), # the index ref genome
                           file.R1,
                           file.R2,
                           "2> /dev/null",
                           "| samtools", "view", "-Sb", 
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



#

# old.file <- list.files( file.path(get.value("align.path"), "Bsaida"), full.names = T)
# 
# new.file <- old.file %>% str_replace(".sorted.bam", ".bam")
# 
# old.file[1:4]
# new.file[1:4]
# 
# file.copy(from = old.file, to = new.file)
# 
# Remove unsorted files that are unecessary 

files.to.remove <- list.files(get.value("align.path"), pattern = ".bam", full.names = T)  %>% 
  str_subset(pattern = "sorted", negate = T)

length(files.to.remove)
files.to.remove[1:20]

for(x in files.to.remove){
  
  file.remove(x)
}

# If you want the big files to be ignore, run the following :

cat("*.bam", "*.tsv", "!.gitignore", sep = "\n",
    file = file.path(get.value("align.path"), ".gitignore"))



# Compute the aligned reads -----------------------------------------------

map.res <- data.frame(ID = character(),
                      total = numeric(),
                      secondary = numeric(),
                      supplementary = numeric(),
                      duplicates = numeric(),
                      mapped = numeric(),
                      mapped_perc = numeric(),
                      paired = numeric(),
                      read1 = numeric(),
                      read2 = numeric(),
                      properly_paired = numeric(),
                      properly_paired_perc = numeric(),
                      twith_itself_mate_mapped = numeric(),
                      singletons = numeric(),
                      singletons_perc = numeric(),
                      twith_mate_mapped_diff_chr = numeric(),
                      twith_mate_mapped_diff_chr_HMQ = numeric(),
                      #Nmappedprim = numeric(),
                      stringsAsFactors = F)

files.to.use <- list.files(get.value("align.path"), pattern = ".stat.tsv", full.names = T)  

for(x in seq_along(files.to.use)){
  
  ID <-   files.to.use[x] %>% str_remove(get.value("align.path")) %>% 
    str_remove(".stat.tsv") %>% 
    str_remove("/")
  
  temp <-  sapply(str_split(readLines(files.to.use[x]), "\t"), `[`,1)  
  map.res[x,] <- c(ID, temp %>% str_remove("%"))
  
  
}  

for(x in 2:ncol(map.res)){
  map.res[,x] <- as.numeric(as.character(map.res[,x]))
  
}  

nrow(map.res)
head(map.res)  
summary(map.res)  

map.res <- map.res %>%  left_join(pop.data, by = c("ID" = "ID_GQ"))

graph1.0 <- map.res %>% #filter(ID != "Dp_2066") %>% 
  ggplot(aes(x = mapped_perc, fill = No_soumission_GQ)) + 
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 95, lty = "dashed") +
  #facet_grid(Espece ~ . , scales = "free")
  labs(x = "Percentage of reads mapped", y = "Density") +
  theme_bw()  

graph1.0

ggsave(filename = file.path(get.value("demulti.log"), "Alignment_density.png"),
       plot = graph1.0,
       width =5, height = 4, units = "in")

graph1.1 <-  map.res %>%  ggplot(aes(y = mapped_perc, x = No_soumission_GQ)) + 
  geom_boxplot() +
  geom_jitter(aes(col = factor(Annee_echantillonnage)),height = 0, alpha = 0.5) +
  #facet_grid(Espece ~ . , scales = "free")
  labs(y = "Percentage of reads mapped") +
  geom_hline(yintercept = 95, lty = "dashed") +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

graph1.1 

ggsave(filename = file.path(get.value("demulti.log"), "Alignment_boxplot.png"),
       plot = graph1.1,
       width =5, height = 4, units = "in")

graph1.2 <- map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=total, y=mapped_perc, col=No_soumission_GQ)) +
  geom_point(alpha = 0.5)+
  geom_hline(yintercept = 95, lty = "dashed") +
  #scale_x_continuous(trans = "log10") + 
  labs(y = "Percentage of reads mapped", x = "N read total") +
  theme_bw() 

graph1.2

graph1.3 <- map.res %>% #filter(Espece_revision %in% c("Pb"))  %>% 
  ggplot(aes(y = mapped_perc, x = No_soumission_GQ)) + 
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = c(95,96), lty = "dashed") +
  geom_jitter(aes(col = total),height = 0, alpha = 0.75) +
  scale_color_distiller(palette = "Spectral", trans = "log10") +
  #facet_wrap(~ Gen_ZONE)
  labs(y = "Percentage of reads mapped") +
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

graph1.3 


graph1.4 <- map.res %>% #filter(Espece_revision %in% c("Pb"))  %>% 
  ggplot(aes(x = mapped, col = No_plaque_envoi)) +
  geom_density() +
  #geom_vline(xintercept = map.res %>% filter(Espece_revision %in% c("Pb")) %>% pull(mapped) %>% median, 
  #           lty = "dashed" )+
  #facet_wrap(~ Gen_ZONE)+
  labs(x = "N reads mapped") +
  theme_bw()

graph1.4

graph1.5 <- map.res %>% #filter(Espece_revision %in% c("Pb"))  %>% 
  ggplot(aes(x = mapped, col = No_soumission_GQ)) +
  geom_density() +
  #geom_vline(xintercept = map.res %>% filter(Espece_revision %in% c("Pb")) %>% pull(mapped) %>% median, 
  #           lty = "dashed" )+
  labs(x = "N reads mapped") +
  #facet_wrap(~ No_soumission_GQ)+
  theme_bw()
graph1.5


## Overall popmap

# Select those that will be used based on alignment
# write.table(map.res %>% filter(mapped_perc >= 96) %>% select(ID) %>% 
#                               mutate(Pop  = "NoPop"), 
#            file = file.path(get.value("info.path"), "popmap_good_align.txt"),
#            quote = FALSE, sep = "\t",
#            row.names = F, col.names = F)

map.res %>% #group_by(Espece_revision) %>% 
  summarise(MedianNread = median(total),
            MedianPerc = median(mapped_perc),
            MinPerc = min(mapped_perc),
            MaxPerc = max(mapped_perc))

plot(map.res[,c("total", "secondary", "mapped", "paired", "properly_paired", "singletons", "twith_itself_mate_mapped", "twith_mate_mapped_diff_chr")])

map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=total, y=secondary, col=No_plaque_envoi)) +
  geom_smooth()+
  #geom_point() +
  #facet_wrap(~Espece_revision)+
  theme_bw()


map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=secondary, y=singletons, col=No_plaque_envoi)) +
  #geom_smooth()+
  geom_point() +
  # facet_wrap(~Espece_revision)+
  theme_bw()


map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=total, y=singletons, col=No_plaque_envoi)) +
  geom_point()

map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=twith_itself_mate_mapped, y=singletons, col=No_plaque_envoi)) +
  geom_point()


 map.res %>% filter(Cat_Sample == "Sample",
                    between(Permapped, 0.95, 0.97)) %>% 
   group_by(Espece, Gen_ZONE) %>% 
   summarise(N = n()) %>% arrange(desc(N)) #%>% pull(N) %>% sum()

hist(map.res$total) 

quantile(map.res$mapped, c(0.01,0.05,0.1))# %>% summary()

# map.res %>%  %>% summarise(Mean = mean(Permapped))

#write_csv(map.res, file = "./02_Results/Mapping_GarMor3_20230817.csv")

# # Stacks - Testing parameters - De-novo ----------------------------------------
# 
# # parameters 
# m.para <- c(3) # 3 is the recommand value .... but can test 4 and 5
# M.para <- c(1:8) # Maybe test the 0
# n.delta <- c(-1, 0, 1)# n = M
# 
# get.value("test.denovo.path")
# get.value("test.denovo.log")
# 
# file.exists(get.value("test.denovo.path"))
# file.exists(get.value("test.denovo.log"))
# file.exists(file.path(get.value("info.path"), "popmap.test_samples.txt"))
# #dir.create(get.value("test.denovo.path"))
# 
# # WILL RUN ONLY IF THE FILE DOESN'T EXIST!!!
# # Check if it's paired or not ...
# 
# for(m in m.para){
#   for(M in M.para){
#     for(nn in n.delta){
#       
#       # last missing parameters
#       n.para <- M + nn
#       
#       # Directory name
#       test.dir <-  file.path(get.value("test.denovo.path"), paste0("stack.m",m,".M",M,".n",n.para))
#       
#       # Run the code only if the directory doesn't exist
#       if(file.exists(test.dir) == F){
#         
#         cat("\nProcessing: ", test.dir, "\n", sep="")
#         
#         dir.create(test.dir)  
#         
#         # Run denovo_map.pl
#         
#         cmd <- paste("--samples", get.value("demulti.path"),
#                      "--popmap", file.path(get.value("info.path"), "popmap.test_samples.txt"),
#                      "-o", test.dir,
#                      "-M", M,
#                      "-n", n.para,
#                      "-m", m, # default = 3 
#                      "--paired",
#                      "-T", 20,
#                      "-X", "\"populations:-r 0.80\""
#         )
#         
#         A <- system2("denovo_map.pl", cmd, stdout=T, stderr=T)
#         
#         # save a log file 
#         
#         cat(file = file.path(get.value("test.denovo.log"), paste0("stack.m",m,".M",M,".n",n.para,".log" )),
#             cmd, "\n\n",
#             A, # what to put in my file
#             append= F, sep = "\n")
#       } # Close the file exist condition
#       
#     } # close nn (aka n)
#   } # Close M
# } # Close m
# 
# 
# # Extract statistics
# 
# nloci.res <- data.frame(m = character(), M = character(), n = character(),
#                         r80 = character(), n_loci = numeric())
# 
# nsnps.res <- data.frame(m = character(), M = character(), n = character(),
#                         r80 = character(), n_snps = numeric(), n_loci = numeric())
# 
# 
# 
# for(x in list.files(get.value("test.denovo.path"), full.names = T)){
#   
#   # Extract parameters
#   info.param <- x %>% str_remove(get.value("test.denovo.path")) %>% 
#     str_remove("/stack.")
#   
#   m.param <- sapply(str_split(info.param, "\\."), `[`)[[1]] %>% str_remove("m")
#   M.param <- sapply(str_split(info.param, "\\."), `[`)[[2]] %>% str_remove("M")
#   n.param <- sapply(str_split(info.param, "\\."), `[`)[[3]] %>% str_remove("n")
#   
#   # N loci
#   
#   cmd <- paste(file.path(x, "populations.log.distribs"),
#                "samples_per_loc_postfilters")
#   
#   res <- system2("stacks-dist-extract", cmd, stdout = T)
#   
#   data <- res[c(-1, -2)] %>% str_split(pattern = "\t")
#   
#   data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
#   names(data) <- c("n_samples", "n_loci")
#   
#   nloci.res <- bind_rows(nloci.res, 
#                          data.frame(m = m.param, M = M.param, n = n.param,
#                                     r80 = "yes", n_loci = as.numeric(as.character(data$n_loci)) %>%  sum() )
#   )
#   
#   # N SNP / locus
#   
#   cmd <- paste(file.path(x, "populations.log.distribs"),
#                "snps_per_loc_postfilters")
#   
#   res <- system2("stacks-dist-extract", cmd, stdout = T)
#   
#   data <- res[c(-1, -2)] %>% str_split(pattern = "\t")
#   
#   data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
#   names(data) <- c("n_snps", "n_loci")
#   
#   data$n_snps <- as.numeric(as.character(data$n_snps))
#   data$n_loci <- as.numeric(as.character(data$n_loci))  
#   
#   data$m   <- m.param
#   data$M   <- M.param
#   data$n   <- n.param
#   data$r80 <- "yes"  
#   
#   nsnps.res <- bind_rows(nsnps.res, data)   
#   
# }
# 
# nloci.res <- nloci.res %>% left_join(nsnps.res %>% filter(n_snps == 0) %>% select(m,n,M, no_div = n_loci), 
#                                      by = c("m", "M", "n")) %>% 
#   mutate(n_loci_poly = n_loci - no_div)
# nloci.res
# nsnps.res
# 
# # Le meilleur paramètre
# nloci.res %>% select(m, M, n, n_loci_poly) %>% arrange(desc(n_loci_poly)) %>% head(10)
# 
# graph1.1 <- nloci.res %>% gather(n_loci, n_loci_poly, key = "loci", value = "n_loci") %>%  
#   mutate(n.rel = ifelse(M==n, "M", ifelse(n>M, "M + 1", "M - 1")),
#          grouped = paste(loci, n.rel)) %>% 
#   ggplot(aes(x = M, y = n_loci, group = grouped, col = n.rel)) +
#   geom_line(aes(lty = loci), show.legend = F) +
#   geom_point() +
#   labs(y = "No. of loci\nshared by 80% of samples", 
#        colour = "n") + 
#   facet_grid(.~m, labeller = label_both) +
#   guides(colour = guide_legend(title.hjust = 0.5))+
#   theme_bw()
# 
# graph1.1
# 
# graph1.2 <- nsnps.res %>% mutate(n_snps_tot = n_snps * n_loci) %>% 
#   group_by(m, M, n) %>% 
#   summarise(n_snps = sum(n_snps_tot)) %>% 
#   mutate(n.rel = ifelse(M==n, "M", ifelse(n>M, "M + 1", "M - 1"))) %>% 
#   ggplot(aes(x = M, y = n_snps, group = n.rel, col = n.rel)) +
#   geom_line() +
#   geom_point() +
#   labs(y = "No. of snps\nshared by 80% of samples", colour = "n") + 
#   facet_grid(.~m, labeller = label_both) +
#   guides(colour = guide_legend(title.hjust = 0.5))+
#   theme_bw()
# 
# graph1.2
# 
# graph1.3 <- nsnps.res %>% mutate(n = ifelse(M==n, "M", ifelse(n>M, "M + 1", "M - 1"))) %>% 
#   filter(m == 3, n %in% c("M + 1", "M", "M - 1")) %>% 
#   group_by(m, M, n) %>% 
#   mutate(n_snps = ifelse(n_snps > 10, 11, n_snps),
#          p_loci = n_loci / sum(n_loci)) %>% 
#   ggplot(aes(x = n_snps, y = p_loci, fill = M)) + 
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(y = "Percentage of loci") + 
#   facet_grid(n ~ m, labeller = label_both)+
#   guides(colour = guide_legend(title.hjust = 0.5))+
#   theme_bw()
# 
# graph1.3
# 
# ggsave(filename = file.path(get.value("test.denovo.log"), "NLoci_test.png"),
#        plot = graph1.1,
#        width = 7.5, height = 4, units = "in")
# 
# ggsave(filename = file.path(get.value("test.denovo.log"), "NSnps_test.png"),
#        plot = graph1.2,
#        width = 7.5, height = 4, units = "in")
# 
# ggsave(filename = file.path(get.value("test.denovo.log"), "NSnpsByLocus_test.png"),
#        plot = graph1.3,
#        width = 7.5, height = 4, units = "in")
# 


# Stacks - Testing parameters - Ref-genome --------------------------------

# This part will run the ref_map.pl pipeline with a subset of samples,
# just to be sure that everything is alright

cmd <- paste("--samples", file.path(get.value("align.path"), "Bsaida"),
             "--popmap", file.path(get.value("info.path"), "popmap.test_samples.txt"),
             "-o", file.path(get.value("test.ref.path"), "Bsaida"),
             "-X", "\"gstacks:-S .sorted.bam\"", 
             "-T", 8,
             "-X", "\"populations:-r 0.80 --vcf --write-single-snp\""
             # en espérant que ça fonctionne ...
)

cmd

A <- system2("ref_map.pl", cmd, stdout=T, stderr=T)
A

# Stats

# Extract statistics

cmd1 <- paste(file.path(get.value("test.ref.path"), "populations.log.distribs"),
              "samples_per_loc_postfilters")

res1 <- system2("stacks-dist-extract", cmd1, stdout = T)

data.int <- res1[c(-1, -2)] %>% str_split(pattern = "\t")

data.int <-  data.frame(matrix(unlist(data.int), nrow= length(data.int), byrow = T))
names(data.int) <- c("n_samples", "n_loci")

n_loci <- as.numeric(as.character(data.int$n_loci)) %>%  sum() 

# N SNP / locus

cmd2 <- paste(file.path(get.value("test.ref.path"), "populations.log.distribs"),
              "snps_per_loc_postfilters")

res2 <- system2("stacks-dist-extract", cmd2, stdout = T)

data.int <- res2[c(-1, -2)] %>% str_split(pattern = "\t")

data.int <-  data.frame(matrix(unlist(data.int), nrow= length(data.int), byrow = T))
names(data.int) <- c("n_snps", "n_loci")

data.int$n_snps <- as.numeric(as.character(data.int$n_snps))
data.int$n_loci <- as.numeric(as.character(data.int$n_loci))  

no_div <- data.int %>% filter(n_snps == 0) %>% pull(n_loci)  %>% as.numeric()
n_loci_poly <- n_loci - no_div

cat("\nThere is", n_loci_poly, "polymorphic loci (r80) out of", n_loci, "loci")


# If you want the big files to be ignore, run the following :

cat("catalog.*", "*.tsv", "!.gitignore", sep = "\n",
    file = file.path(get.value("test.ref.path"), ".gitignore"))



vcf.path <- ("00_Data/04b_Test.ref/populations.snps.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

library(remotes)
#remotes::install_github("biodray/QuickPop")
library(QuickPop)
library(adegenet)
gl.data  <- vcfR::vcfR2genlight(vcf.data) 
pop(gl.data) <- data.frame(ID_GQ = indNames(gl.data)) %>% 
  left_join(pop.data) %>% pull(Espece)

table(pop(gl.data))

pca.test  <- glPca(gl.data, center = TRUE, scale = FALSE,  
                   parallel = TRUE, n.core = 20, nf = 1000)

gPCA <- pca.test %>% QuickPop::pca_scoretable(naxe = 3) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col =Region_echantillonnage)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data)), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.test)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.test)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA

ggsave(filename = file.path(get.value("test.ref.log"), "PCA_test_StackRef.png"),
       plot = gPCA,
       width = 5, height = 4, units = "in")

# If you want the big files to be ignore, run the following :
cat("catalog.*", "*.tsv", "*.vcf", "!.gitignore", sep = "\n",
    file = file.path(get.value("test.ref.path"), ".gitignore"))

# # Stacks - Testing parameters - Comparing de-novo VS ref ------------------
# 
# # N SNP / locus - REF m3 M2 n3
# 
# cmd2 <- paste(file.path(get.value("test.denovo.path"), "stack.m3.M2.n3", "gstacks.log.distribs"),
#               "effective_coverages_per_sample")
# 
# res2 <- system2("stacks-dist-extract", cmd2, stdout = T)
# 
# data2 <- res2[c(-1, -2, -3)] %>% str_split(pattern = "\t")
# 
# data2 <-  data.frame(matrix(unlist(data2), nrow= length(data2), byrow = T))
# names(data2) <- c("sample", "n_loci", "n_used_fw_reads", "mean_cov", "mean_cov_ns")
# data2 <- data2 %>% mutate(Method = "denovo")
# 
# # N SNP / locus - denovo
# 
# cmd1 <- paste(file.path(get.value("test.ref.path"), "gstacks.log.distribs"),
#               "effective_coverages_per_sample")
# 
# res1 <- system2("stacks-dist-extract", cmd1, stdout = T)
# 
# data1 <- res1[c(-1, -2, -3)] %>% str_split(pattern = "\t")
# 
# data1 <-  data.frame(matrix(unlist(data1), nrow= length(data1), byrow = T))
# names(data1) <- c("sample", "n_loci", "n_used_fw_reads", "mean_cov", "mean_cov_ns")
# data1 <- data1 %>% mutate(Method = "ref")
# 
# data.test <- rbind(data1, data2)
# 
# data.test 
# 
# graph1.4 <- data.test  %>% select(sample, mean_cov_ns, Method) %>% 
#   spread("Method", mean_cov_ns)  %>% 
#   left_join(data, by = c("sample" = "ID_GQ")) %>%
#   ggplot(aes(x = as.numeric(as.character(ref)), y = as.numeric(as.character(denovo)), col = Espece)) +
#   geom_point() +
#   geom_abline(slope = 1, lty = "dashed", col = "darkgray") +
#   #facet_wrap(~POP, nrow = 2) +
#   labs(x = "Mean coverage - reference genome", y = "Mean coverage - denovo") +
#   scale_x_continuous(limits = c(3,20)) +
#   scale_y_continuous(limits = c(3,20)) +
#   theme_bw() +
#   theme(legend.position = "right")
# 
# graph1.4
# 
# ggsave(filename = file.path(get.value("stacks.ref.log"), "CoverageBySample_ComparisonRefvsDenovo.png"),
#        plot = graph1.4,
#        width = 4, height = 4, units = "in")
# 
# 
# 
# 
# 


# Gstakcs - Ref - discovered SNPs -----------------------------------------

cmd <- paste("-I", get.value("align.path"),
             "-M", file.path(get.value("info.path"), "popmap_good_align.txt"),
             "-O", get.value("stacks.ref.path"),
             "-S", ".sorted.bam",
             "-t", 8)

A <- system2("gstacks", cmd, stdout=T, stderr=T)

cat(file = file.path(get.value("stacks.ref.log"), "gstacks.ref.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")



# Check the distribution of ...


# N SNP / locus

cmd <- paste(file.path(get.value("stacks.ref.path"), "gstacks.log.distribs"),
             "effective_coverages_per_sample")

res <- system2("stacks-dist-extract", cmd, stdout = T)

data <- res[c(-1, -2, -3)] %>% str_split(pattern = "\t")

data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
names(data) <- c("sample", "n_loci", "n_used_fw_reads", "mean_cov", "mean_cov_ns")

library(ggridges)

#write_csv(data, file.path(get.value("stacks.ref.log"), "AllIndividuals_NreadsNloci.csv"))


graph2.1 <- data %>% left_join(pop.data, by = c("sample" = "ID_GQ")) %>% 
  filter(as.numeric(as.character(mean_cov_ns)) < 100) %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), y  = No_plaque_envoi, fill = ..x..)) +
  ggridges::geom_density_ridges_gradient(alpha = 1/2, scale = 3) +
  scale_fill_viridis_c(name = "Mean coverage") +
  
  geom_vline(xintercept = c(5,10), lty = "dashed", col = "black") +
  scale_x_continuous(breaks = c(0,5,10,20,40,60,80, 100))  + #facet_wrap(~ Objective) +
  labs(x = "Mean coverage", title = "Mean coverage post gstacks in Mackerels") + 
  # theme_ridges() +
  theme_bw() +
  theme(legend.position = "right")

graph2.1

ggsave(filename = file.path(get.value("stacks.ref.log"), "CoverageByPlate_Aug2023.png"),
       plot = graph2.1,
       width = 7.5, height = 4, units = "in")

graph2.2 <- data %>% left_join(pop.data, by = c("sample" = "ID_GQ")) %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), y = as.numeric(as.character(n_loci)), col = No_plaque_envoi)) +
  geom_point(alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  geom_vline(xintercept = c(5,10), lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = 50000, lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(x = "Mean coverage", y = "N loci") + 
  theme_bw() +
  theme(legend.position = "right")

graph2.2

ggsave(filename = file.path(get.value("stacks.ref.log"), "CoverageVSloci_May2023.png"),
       plot = graph2.2,
       width = 7, height = 4, units = "in")

graph2.3 <- data %>% left_join(pop.data, by = c("sample" = "ID_GQ")) %>% 
  ggplot(aes(y = as.numeric(as.character(mean_cov_ns)), x = Region_echantillonnage, col = Region_echantillonnage)) +
  geom_violin(col = "black") +  
  geom_jitter(height = 0, alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  #geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = c(5,10), lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(y = "Mean coverage", y = "") + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

graph2.3

ggsave(filename = file.path(get.value("stacks.ref.log"), "Coverage_Aug2023.png"),
       plot = graph2.3,
       width = 7, height = 4, units = "in")





