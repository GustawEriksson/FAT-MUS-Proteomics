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

Enrichment.Analysis <- function(DEP.table, selected.db, output.dir, name) {
  
  setEnrichrSite("Enrichr")
  DEP.signif = DEP.table$name[DEP.table$significant == TRUE]
  enriched <- enrichr(DEP.signif, selected.db)
  
  for (i in names(enriched)) {
    print(i)
    Enrich.plot = plotEnrich(enriched[[i]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",
                             title = i)
    ggsave(paste0(output.dir, name, "_",i, "_enrichR.pdf"), plot = Enrich.plot)
    writexl::write_xlsx(x = enriched[[i]], path = paste0(output.dir, name, "_",i, "_enrichR_table.xlsx"))
    
  }
}