# Author: Gustaw Eriksson
# Date: 2021-09-16
# Description: Quick and dirty analysis of proteomics data of 
# fat and muscle received from Anna Benrick

# CHECK ADJUSTED P-VALUE IN DEP OUTPUT. SOME P.ADJ IS LOWER THAN
# NON-ADJUSTED P-VALUE.

#library(DEqMS)
library(matrixStats)
library(ggplot2)
library(cowplot)
library(readxl)
#library(ComBat)
#library(gmm)
library(DEP)
library(dplyr)
#library(ggpubr)
library(ggfortify)
library(reshape2)
#library(EnhancedVolcano)
#library(mice)
library(sva)
library(data.table)
library(enrichR)

# Set work directory to directory where script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load function scripts in /bin
source(file = "bin/DEP-Preprocessing.R")
source(file = "bin/DEP-ComBat.Batch.Correction.R")
source(file = "bin/DEP_Differential.Expression.Analysis.R")
source(file = "bin/DEP.Enrichment_Analysis.R")


# Setting parameters

  # Is VSN data normalisation to be done
VSN_norm = FALSE
  # Set alpha threshold (p-value significans cut-off)
alpha_set = 0.05
  # Set log2 fold change threshold 
log2fc_set = log2(1.5)

# Load either fat or muscle data and reformatting the data
Load_muscle = TRUE
Load_fat = FALSE

# Set output directory in project fold
data_output = "Output/Muscle_No-VSN_TopValue_220217/"

if (Load_muscle == TRUE) {
  data_raw = readxl::read_xlsx("Data/MUSCLE_TMT_Intensity_edited.xlsx")
  #data_raw = subset(data_raw, select = -c(3,18:47, 63, 79))
  #data_raw_orig = cbind(data_raw[,1:16], data_raw[,17:46][order(sub('.*m ', '', colnames(data_raw[,17:46])))])
  #colnames(data_raw)[17:46] = c(paste0("TMT_m", seq(30)))
  
  data_raw = subset(data_raw, select = -c(3,18:47, 63, 79))
  data_raw = cbind(data_raw[,1:16], data_raw[,17:46][order(sub('.*m ', '', colnames(data_raw[,17:46])))])
  colnames(data_raw)[17:46] = c(paste0("TMT_m", seq(30)))
  
  # Setting muscle parameters
  exp.design.matrix = "Data/Muscle_experimental_design.xlsx"
  tissue.type = "Muscle"
  
} else if (Load_fat == TRUE) {
  data_raw = readxl::read_xlsx("Data/FAT_TMT_Intensity_edited.xlsx")
  data_raw = subset(data_raw, select = -c(3,18:47, 63, 79))
  data_raw_orig = cbind(data_raw[,1:16], data_raw[,17:46][order(sub('.*f ', '', colnames(data_raw[,17:46])))])
  colnames(Muscle_raw)[17:46] = c(paste0("TMT_f", seq(30)))
  
  # Setting fat parameters
  exp.design.matrix = "Data/Fat_experimental_design.xlsx"
  tissue.type = "Fat"
}

# Preprocessing fo the data to rename duplicates, filter and impute 
# missing values, and plot the results
data.preprocess = DEP_Preprocess(P.matrix = data_raw, 
                                   exp.design = exp.design.matrix, 
                                   output.dir = data_output, 
                                   name = tissue.type,
                                   vsn = VSN_norm, 
                                   keep.duplicates = TRUE)

# Batch correct the data using ComBat batch correction
data.corr = Batch.Correction(P.matrix = data.preprocess$name, 
                            batch.order = c(rep(1,5), rep(2,5), rep(1,5), rep(2,5), rep(1,5), rep(2,5)), 
                            batch.colour = c(rep("Batch 1", 5), rep("Batch 2", 5), rep("Batch 1", 5),
                                             rep("Batch 2", 5), rep("Batch 1", 5), rep("Batch 2", 5)), 
                            batch.shape = c(rep("Control", 10), rep("PCOS w0", 10), rep("PCOS w5", 10)), 
                            name = tissue.type,
                            output.dir = data_output)

# Do differential expression analysis of protein data
data.DEP = DEP_analysis(P.corr = data.corr,
                       ids = data.preprocess$ID,
                       name = tissue.type,
                       output.dir = data_output,
                       exp.design = exp.design.matrix,
                       alpha = alpha_set,
                       log2fc = log2fc_set)

#data.table_generator = DEP_table.generator(output.dep = data.DEP,
#                                           output.dir = data_output,
#                                           name = tissue.type)

data.w5_w0 = Generate_DEP_table(x = data.DEP$PCOS_w5_w0_p.adj, name = paste0(tissue.type, "_PCOS_w0_w5"),
                               output.dir = data_output)
data.Ctrl_w0 = Generate_DEP_table(x = data.DEP$`PCOS_w0-Ctrl_p.adj`, name = paste0(tissue.type, "_Ctrl_PCOS-w0"),
                                 output.dir = data_output)

# Enrichment analysis with enrichR
#setEnrichrSite("Enrichr") # Human genes
#dbs <- listEnrichrDbs()
#dbs_selected = c("KEGG_2021_Human", "GO_Biological_Process_2021", "GO_Cellular_Component_2021",
#                 "GO_Molecular_Function_2021")

#DEP.table = Fat_w5_w0
#selected.db = dbs_selected
#output.dir = Fat_output
#name = "Fat_PCOS_w0_w5"

data.w5_w0 = Enrichment.Analysis(DEP.table = data.w5_w0, 
                                   name = paste0(tissue.type, "_PCOS_w0_w5"),
                                   output.dir = data_output, 
                                   selected.db = c("KEGG_2021_Human", "GO_Biological_Process_2021", "GO_Cellular_Component_2021",
                                                   "GO_Molecular_Function_2021"))

data.Ctrl_w0 = Enrichment.Analysis(DEP.table = data.Ctrl_w0, 
                                     name = paste0(tissue.type, "_Ctrl_PCOS-w0"),
                                     output.dir = data_output, 
                                     selected.db = c("KEGG_2021_Human", "GO_Biological_Process_2021", "GO_Cellular_Component_2021",
                                                     "GO_Molecular_Function_2021"))

# Output and save batch corrected table
data.corr_df = setDT(data.corr, keep.rownames = TRUE)
colnames(data.corr_df)[1] = "Protein" 
writexl::write_xlsx(x = data.corr_df, 
                    path = paste0(data_output, tissue.type,"_ComBat_Protein_table.xlsx"))

# Move plots to sharing directory
#Share.dir <- function(output.dir) {
#  dir.create(path = paste0(output.dir, "Share"), showWarnings = FALSE)
#  copy.list = list.files(path = output.dir, full.names = TRUE)[c(3:11, 27:35, 40:43)]
#  file.copy(from = copy.list, to = paste0(output.dir, "Share"), overwrite = TRUE)
#}

#data.share = Share.dir(output.dir = data_output)

