#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()

BiocManager::version()

#BiocManager::install("vcfR")

library(vcfR)

#BiocManager::install(c("DESeq2", "recount3"))

library(DESeq2)
library(recount3)
