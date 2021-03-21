if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!requireNamespace("devtools")) BiocManager::install(c("devtools"))
if(!requireNamespace("robustbase"))BiocManager::install(c("robustbase"))
library(devtools)
library(robustbase)
devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
if(!requireNamespace("SummarizedExperiment"))BiocManager::install(c("SummarizedExperiment"))
if(!requireNamespace("maftools"))BiocManager::install(c("maftools"))
if(!requireNamespace("arsenal"))install.packages(c("arsenal"))
if(!requireNamespace("survival"))install.packages(c("survival"))
if(!requireNamespace("survminer"))install.packages(c("survminer"))

library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
library(arsenal)
library(survival)
library(survminer)

mutation <- GDCquery_Maf(tumor = "BRCA",save.csv=TRUE, pipeline="mutect2")
maf_dataframe = read.maf(mutation)

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query) #only need this line of code ONCE to download the data
sum_exp <- GDCprepare(query)


patient_data <- colData(sum_exp)
patient_ages <- patient_data$paper_age_at_initial_pathologic_diagnosis
#here we are creating a NEW COLUMN in patient_data called "age_category"
#NOTE: This will NOT be added to colData(sum_exp). Instead it will only be added to the patient_data data table.
patient_data$age_category = ifelse(patient_ages < 40, "Young", ifelse(patient_ages >= 60, "Old", "Mid"))
#The ifelse() form is: ifelse( condition, action when condition is true, action when condition is false ). Here we have two ifelse() embedded together

patient_data <- patient_data[, c("barcode", "age_category")]

patient_data

#get shortened patient barcodes so that we can compare with
short_maf <- substr(maf_dataframe@clinical.data$Tumor_Sample_Barcode, 1,12)

#create a new column in maf_dataframe
maf_dataframe@clinical.data$short_barcodes <- short_maf

#extract age category information for each barcode that we have
maf_ages <- patient_data[short_maf, "age_category"]

maf_dataframe@clinical.data$Ages <- maf_ages

#Extract codes for each age group
young_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Young",]

mid_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Mid",]

old_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Old",]

#create maf subsets for each age group
young_maf <- subsetMaf(maf_dataframe, tsb = young_codes$Tumor_Sample_Barcode)

mid_maf <- subsetMaf(maf_dataframe, tsb = mid_codes$Tumor_Sample_Barcode)

old_maf <- subsetMaf(maf_dataframe, tsb = old_codes$Tumor_Sample_Barcode)

#get oncoplots
oncoplot(maf = maf_dataframe, draw_titv = TRUE, top = 3)
oncoplot(maf = young_maf, draw_titv = TRUE, top = 3)
oncoplot(maf = mid_maf, draw_titv = TRUE, top = 3)
oncoplot(maf = old_maf, draw_titv = TRUE, top = 3)
