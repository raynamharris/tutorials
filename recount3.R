# From https://bioconductor.org/packages/release/bioc/vignettes/recount3/inst/doc/recount3-quickstart.html

library(tidyverse)
library(recount3)
library(cowplot)
library(scales)

# get data
human_projects <- available_projects()
head(human_projects)

subset(human_projects, 
       file_source == "gtex" & 
         project_type == "data_sources")

gtex_heart <- subset(human_projects,
  project == "HEART" & 
    file_source == "gtex" & 
    project_type == "data_sources")

rse_gtex_heart <- create_rse(gtex_heart)
rse_gtex_heart

mycounts <- assays(rse_gtex_heart)$raw_counts %>% as.data.frame()
head(mycounts)
names(mycounts)

mycolData <- colData(rse_gtex_heart) %>% as.data.frame()
head(mycolData)

# check that rows and samples match
rownames(mycolData) == colnames(mycounts)

