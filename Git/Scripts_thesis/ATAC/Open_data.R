#upload all files 
library(openxlsx)
library(readxl)
library(jsonlite)
#BiocManager::install("rtracklayer")
library(rtracklayer)
 


#metadata 
samplesheet <- read.csv("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/atac/atac/data/minkang/sample_input_eralin.csv")

#peak files 

peak_file <- read.table("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/atac/atac/data/peaks/CD8_Tcell_s1_1_202103734_peaks.broadPeak.noblacklist", header = FALSE, sep = "/t")
head(peak_file)


#QC file 
qc <- fromJSON(gzfile("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/atac/atac/data/qc/data/CD8_Tcell_s1_1_202103734.json.gz"))
head(qc)


#signal 

signal <- import("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/atac/atac/data/signal/CD8_Tcell_s1_1_202103734.bigwig")
head(signal)

#QC and meta 
metadata <- read_excel("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/atac/atac/data/QC of Samples for ATACseq (1).xlsx")
head(metadata)

#delete last 3 rows 

#protein

protein <- read.csv("R:/TCR/Biotech/Genmab/Felix Behr/Proteomic Data/Project P1491/Proteome Discoverer Files/P1491-2 T cell stimulation proteomics Felix.pdStudy")
