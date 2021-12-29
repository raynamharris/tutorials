library(tidyverse)
library(ggplot2)
library(DESeq2)

## get data

## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88741

## https://maayanlab.cloud/biojupies/analyze/example

countData <- read.delim("gse88741/GSE88741-expression.txt")
head(countData)
str(countData)

colData <- read.delim("gse88741/GSE88741-metadata.txt")
head(colData)
str(colData)

## make joined df for easy ggplot 2 plotting

df <- countData %>%
  pivot_longer(-gene_symbol, 
               names_to = "Sample_geo_accession",
               values_to = "Count") %>%
  full_join(., colData)
head(df)


## total gene counts 

totals <- df %>%
  group_by(Sample_geo_accession, Sample.Title, Stage, cell.type) %>%
  summarise(Reads = sum(Count))
totals

totals %>%
  ggplot(aes(x = reorder(Sample_geo_accession, Reads),
             y = Reads)) +
  geom_bar(stat = "identity", fill = "#2877b2") +
  scale_y_continuous(labels = scales::label_number_si()) +
  coord_flip() +
  labs(x = "Sample")

## prep for DESeq2 and PCA

rownames(colData) <- colData$Sample_geo_accession
rownames(countData) <- countData$gene_symbol
countData$gene_symbol <- NULL


## for pca
df_pca <- as.data.frame(t(countData)) %>%
  prcomp()
df_pca <- as.data.frame(df_pca$x)
df_pca <- df_pca %>%
  mutate(Sample_geo_accession = row.names(.)) %>%
  full_join(., colData)
head(df_pca)

df_pca %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Stage))

df_pca %>%
  ggplot(aes(x = PC3, y = PC4)) +
  geom_point(aes(color = Stage))

