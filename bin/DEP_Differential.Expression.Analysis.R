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

DEP_analysis <- function(P.corr, ids, output.dir, name, exp.design, alpha = 0.05, log2fc = log2(1.5)) {
  
  P.output = list()
  
  P.corr$name = rownames(P.corr)
  P.corr$ID = rownames(ids)
  P.table = P.corr[,c(32,31, 1:30)]
  P.output[["P.table"]] = P.table
  
  P.corr[,1:30] = 2^P.corr[,1:30]
  P.unique = make_unique(proteins = P.corr, names = "name", ids = "ID")
  P.design <- readxl::read_xlsx(exp.design)
  P.design$label = colnames(P.corr[,1:30])
  P.se = make_se(proteins_unique = P.unique, columns = c(1:30), expdesign = P.design)
  P.output[["P.se.corr"]] = P.se
  
  # PCOS w5 vs. PCOS w0, adjusted p-value
  P.diff <- test_diff(P.se, type = "manual", test = "PCOS_w5_vs_PCOS_w0")
  dep <- add_rejections(P.diff, alpha = alpha, lfc = log2fc)
  P.res = get_results(dep)
  P.output[["PCOS_w5_w0_p.adj"]] = P.res
  plot_volcano(dep, contrast = "PCOS_w5_vs_PCOS_w0", label_size = 2, add_names = TRUE)
  ggsave2(paste0(output.dir,name,"_Volcano_PCOS_w5-w0_p.adj.pdf"))
  plot_pca(dep, x = 1, y = 2, n = 100, point_size = 4, indicate = c("replicate", "condition"))
  ggsave2(paste0(output.dir,name,"_PCA_Labelled_PCOS_w5-w0.pdf"))
  plot_pca(dep, x = 1, y = 2, n = 100, point_size = 4, indicate = "condition")
  ggsave2(paste0(output.dir,name,"_PCA_PCOS_w5-w0.pdf"))
  plot_cor(dep, significant = FALSE, lower = 0, upper = 1, pal = "Reds")
  ggsave2(paste0(output.dir,name,"_PCOS_w5-w0_Pearson_correlation_matrix.pdf"))
  
  
  # PCOS w5 vs. PCOS w0, normal p-value
  P.diff <- test_diff(P.se, type = "manual", test = "PCOS_w5_vs_PCOS_w0")
  P.diff@elementMetadata$PCOS_w5_vs_PCOS_w0_p.adj = P.diff@elementMetadata$PCOS_w5_vs_PCOS_w0_p.val
  
  dep <- add_rejections(P.diff, alpha = alpha, lfc = log2fc)
  P.res = get_results(dep)
  P.output[["PCOS_w5_w0_p.value"]] = P.res
  plot_volcano(dep, contrast = "PCOS_w5_vs_PCOS_w0", label_size = 2, add_names = TRUE)
  ggsave2(paste0(output.dir,name,"_Volcano_PCOS_w5-w0_p.value.pdf"))
  
  # PCOS w0 vs. Control, adjusted p-value
  P.diff <- test_diff(P.se, type = "manual", test = "PCOS_w0_vs_Ctrl")
  dep <- add_rejections(P.diff, alpha = alpha, lfc = log2fc)
  P.res = get_results(dep)
  P.output[["PCOS_w0-Ctrl_p.adj"]] = P.res
  
  plot_volcano(dep, contrast = "PCOS_w0_vs_Ctrl", label_size = 2, add_names = TRUE)
  ggsave2(paste0(output.dir,name,"_Volcano_PCOS_w0-Ctrl_p.adj.pdf"))
  plot_pca(dep, x = 1, y = 2, n = 100, point_size = 4, indicate = c("replicate", "condition"))
  ggsave2(paste0(output.dir,name,"_PCA_Labelled_PCOS_w0-Ctrl.pdf"))
  plot_pca(dep, x = 1, y = 2, n = 100, point_size = 4, indicate = "condition")
  ggsave2(paste0(output.dir,name,"_PCA_PCOS_w0-Ctrl.pdf"))
  plot_cor(dep, significant = FALSE, lower = 0, upper = 1, pal = "Reds")
  ggsave2(paste0(output.dir,name,"_PCOS_w0-Ctrl_Pearson_correlation_matrix.pdf"))
  
  # PCOS w0 vs. Control, normal p-value
  P.diff <- test_diff(P.se, type = "manual", test = "PCOS_w0_vs_Ctrl")
  P.diff@elementMetadata$PCOS_w0_vs_Ctrl_p.adj = P.diff@elementMetadata$PCOS_w0_vs_Ctrl_p.val
  
  dep <- add_rejections(P.diff, alpha = alpha, lfc = log2fc)
  P.res = get_results(dep)
  P.output[["PCOS_w0-Ctrl_p.value"]] = P.res
  plot_volcano(dep, contrast = "PCOS_w0_vs_Ctrl", label_size = 2, add_names = TRUE)
  ggsave2(paste0(output.dir,name,"_Volcano_PCOS_w0-Ctrl_p.value.pdf"))
  
  return(P.output)
}
########
# This function is not used
DEP_table.generator <- function(output.dep, output.dir, name) {
  pdf(file = paste0(output.dir, name, "_meanSDplot_Post-Corr.pdf"))
  meanSdPlot(output.dep$P.se.corr)
  dev.off()
  plot_normalization(output.dep$P.se.corr)
  ggsave(paste0(output.dir, name,"_boxplot_transformation_post.corr.pdf"))
  plot_imputation(output.dep$P.se.corr)
  ggsave(paste0(output.dir, name,"_Imputation_plot_post.corr.pdf"))
  writexl::write_xlsx(x = output.dep$`PCOS_w0-Ctrl_p.adj`, 
                      path = paste0(output.dir, name,"_DEP_PCOS-w0_vs._Ctrl.xlsx"))
  writexl::write_xlsx(x = output.dep$PCOS_w5_w0_p.adj, 
                      path = paste0(output.dir, name,"_DEP_PCOS-w5_vs._w0.xlsx"))
}
############

Generate_DEP_table <- function(x, name, output.dir, alpha = alpha_set, log2fc = log2fc_set) {
  
  x.out = x[,c(1:4, 7)]
  colnames(x.out) = c("name", "ID", "p.value", "p.adj", "Log2FC")
  
  x.out$significant = (x.out$p.value < alpha & x.out$Log2FC < -log2fc) | 
    (x.out$p.value < alpha & x.out$Log2FC > log2fc)
  
  writexl::write_xlsx(x = x.out, path = paste0(output.dir, name, "_DEP_table.xlsx"))
  
  return(x.out)
  
}
