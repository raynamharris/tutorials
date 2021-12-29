# From https://bioconductor.org/packages/release/bioc/vignettes/recount3/inst/doc/recount3-quickstart.html

library(tidyverse)
library(recount3)
library(biomaRt)
library(Rtsne)
library(cowplot)


# get GTEx heart data
human_projects <- available_projects(organism = "human")
head(human_projects)

subset(human_projects, 
       file_source == "gtex" & 
         project_type == "data_sources")

gtex_heart <- subset(human_projects,
                     project == "HEART"  & 
                       file_source == "gtex" & 
                       project_type == "data_sources" )

rse_gtex_heart <- create_rse(gtex_heart)
rse_gtex_heart



# format data for DESEq2
countData <- assays(rse_gtex_heart)$raw_counts %>% as.data.frame()
colData <- colData(rse_gtex_heart) %>% as.data.frame()

# check that rows and samples match
rownames(colData) == colnames(countData)

# variables
dim(colData)
dim(countData)

head(rownames(colData))
(colnames(colData))
head(rownames(countData))
head(colnames(countData))


# subset data
# include only those with SRA accession numbers
# exclude miRNA samples

colData <-  colData %>%
  filter(gtex.run_acc != "NA",
         gtex.smnabtcht != "RNA isolation_PAXgene Tissue miRNA") %>%
  dplyr::select(external_id, study, gtex.run_acc, 
                gtex.age, gtex.smtsd)
head(colData)

ggplot(colData, aes(x = gtex.age, fill = gtex.smtsd)) +
  geom_bar(stat = "count", position = "dodge") +
  labs(#x = "Age", y = "Count", fill = "Tissue",
       subtitle = "GTEx data obtained using recount3 ")

# get countdata for this subset of colData

## colData and countData must contain the exact same samples. 
savecols <- as.character(rownames(colData)) #select the rowsname 
savecols <- as.vector(savecols) # make it a vector
countData <- countData %>% dplyr::select(one_of(savecols)) # select just the columns 
head(countData)[1:5]  


# check that rows and samples match
rownames(colData) == colnames(countData)

## get and clean ensemble ids

ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_info <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),  
                   mart = ensembl) %>% 
  mutate_all(na_if, "") %>%
  drop_na(.) %>%
  unique(.) %>%
  mutate(ensembl_gene_id = paste0(ensembl_gene_id, ".1", sep = ""))

head(gene_info)
tail(gene_info)


## make long for easier ggplotting and join with gene info

countData_long <- countData %>%
  mutate(ensembl_gene_id = rownames(.)) %>%
  pivot_longer(-ensembl_gene_id, 
               names_to = "external_id", values_to = "counts") %>%
  inner_join(gene_info, .,by = "ensembl_gene_id") %>% 
  full_join(colData, ., by = "external_id") %>%
  arrange(desc(counts))
head(countData_long)
tail(countData_long)

# plot your favorite gene

a <- countData_long %>%
  filter( hgnc_symbol == "MT-CO2") %>%
  ggplot(aes(x = gtex.age, y = counts, 
             fill = gtex.smtsd)) +
  geom_boxplot() +
  geom_point()  +
  facet_wrap(~gtex.smtsd) +
  scale_y_log10() +
  labs(y = 'MT-CO2 counts', x = "Age") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1))
a


# widen for tsne

countData_long_wide <- countData_long %>%
  dplyr::select(-hgnc_symbol) %>%
  pivot_wider(id_cols = external_id:gtex.smtsd,
              names_from = ensembl_gene_id, 
              values_from = counts,
              values_fn = sum)
head(countData_long_wide)
colnames(countData_long_wide)

tsne_data <- countData_long_wide[ ,6:306] # We are sub-setting IR object such as to include 'all rows' and columns 1 to 4.
tsne_samples <- countData_long_wide[ ,1:5] # We are sub-setting IR object such as to include 'all rows' and column 5.


## Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_results <- Rtsne(tsne_data, perplexity=30, 
                      check_duplicates = FALSE) 

tsne_results_samples <- as.data.frame(tsne_results$Y) %>%
  cbind(tsne_samples, .)
head(tsne_results_samples) 

b <- tsne_results_samples %>%
  ggplot(aes(x = V1, y = V2, color = gtex.smtsd)) +
  geom_point() +
  theme(legend.position = c(.7,.9),
        legend.direction = "vertical") +
  labs(x = "tSNE dimention 1", 
       y = "tSNE dimention 2", 
       color = "GTEx Tissue")

# plot fav gene and tsne
p <- plot_grid(a,b)
p 

png("../images/recount3-gtex.png")
print(p)
dev.off()

