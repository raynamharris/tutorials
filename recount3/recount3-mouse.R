# From https://bioconductor.org/packages/release/bioc/vignettes/recount3/inst/doc/recount3-quickstart.html
# Manuscript https://elifesciences.org/articles/14997

library(tidyverse)
library(recount3)
library(biomaRt)
library(Rtsne)
library(cowplot)



## recount 3

mouse_projects <- available_projects(organism = "mouse")
head(mouse_projects)

SRP066161 <- subset(mouse_projects, project == "SRP066161" )

myrse <- create_rse(SRP066161)
myrse

## format data for DESEq2

countData <- assays(myrse)$raw_counts %>% as.data.frame()
colData <- colData(myrse) %>% as.data.frame()

# check that rows and samples match
rownames(colData) == colnames(countData)

# variables
dim(colData)
dim(countData)

head(rownames(colData))
head(rownames(countData))
head(colnames(countData))


# subset data
# include only those with SRA accession numbers
# exclude miRNA samples

colData <-  colData %>%
  dplyr::select(external_id:sra.submission_acc, sra.experiment_title,
                sra.run_total_bases) %>%
  separate(sra.experiment_title, into = c("GSM", "sample"), sep = " ") %>%
  separate(sample, into = c("tissue", "rep"), sep = "_")
head(colData)


mylevels <- c("dgd", "dgv", "ca4", "ca3d", "ca3v", "ca2", "ca1d", "ca1v")
mylables <- c("DG dors.", "DG vent.", "CA4", "CA3 dors.", "CA3 vent.", 
              "CA2", "CA1 dors.", "CA1 vent.")
mycolors <- c("#d7322d", "#d7322d", "#612150", "#3b841e", "#3b841e", 
              "#3d98b2", "#0a1550", "#0a1550")
myalpha <- c(1, 0.7, 1, 1, 0.7, 1, 1, 0.7)
myalpha2 <- c(1, 0.5, 1, 1, 0.5, 1, 1, 0.5)


ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

gene_info <- getBM(attributes=c('ensembl_gene_id','mgi_symbol'),  
                   mart = ensembl) %>% 
  mutate_all(na_if, "") %>%
  drop_na(.) %>%
  unique(.) %>%
  separate(ensembl_gene_id, into = c("ensembl_gene_id", NA), sep = "\\.")
head(gene_info)

countData_long <- countData %>%
  mutate(ensembl_gene_id = rownames(.)) %>%
  separate(ensembl_gene_id, into = c("ensembl_gene_id", NA), sep = "\\.") %>%
  pivot_longer(-ensembl_gene_id, 
               names_to = "external_id", values_to = "counts") %>%
  right_join(gene_info, .,by = "ensembl_gene_id") %>% 
  full_join(colData, ., by = "external_id") %>%
  arrange(desc(counts)) %>%
  mutate(tissue = factor(tissue, levels = mylevels))
head(countData_long)
tail(countData_long)

# screen shots to recreate from https://elifesciences.org/articles/14997
fig_broad <- png::readPNG("images/recount3-broad.png")
fig_broad <- ggdraw() +  draw_image(fig_broad, scale = 1)

fig_specific <- png::readPNG("images/recount3-specific.png")
fig_specific <- ggdraw() +  draw_image(fig_specific, scale = 1)




fig2b <- function(mygene, mytitle, myylab){
  p <- countData_long %>%
    filter( mgi_symbol  == mygene) %>%
    group_by(tissue) %>%
    summarize(meancounts = mean(counts),
              sdcounts = sd(counts)) %>%
    ggplot(aes(x = tissue, y = meancounts, 
               fill = tissue, alpha = tissue  )) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin=meancounts-sdcounts, 
                      ymax=meancounts+sdcounts), width=.2) +
    scale_y_continuous(labels = scales::label_number_si(accuracy = 0.1)) +
    scale_alpha_manual(values = myalpha) +
    labs(y = myylab, x = "Cells",
         subtitle = mygene,
         title = mytitle) +
    theme_linedraw(base_size = 15) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, 
                                     hjust = 1),
          plot.subtitle = element_text(face = "italic")) +
    scale_fill_manual(values = mycolors) +
    scale_x_discrete(labels = mylables)
  return(p)
}

a <- fig2b("Prox1", "Granule cells", "\nRNA-Seq")
b <- fig2b("Dkk3", "Non-granule cells", " ")
c <- fig2b("Calb2", "Mossy cells", " ")
d <- fig2b("Ociad2", "All pyramids", " ")
e <- fig2b("Cacng5", "CA2 pyramids", " ")
f <- fig2b("Fibcd1", "CA1 pyramids", " ")

B1 <- plot_grid(a,b,c, d, e, f, nrow =1 ,
                rel_widths = c(1.1,1,1,1,1,1))


g <- fig2b("Pdzd2", "Granule cells, dors.", "\nRNA-Seq")
h <- fig2b("Tox3", "Granule cells, vent.", " ")
i <- fig2b("Iyd", "CA3 pyramids, dors.", " ")
j <- fig2b("Coch", "CA3 pyramids, vent.", " ")
k <- fig2b("Wfs1", "CA1 pyramids, dors.", " ")
l <- fig2b("Dcn", "CA1 pyramids, vent.", " ")

B2 <- plot_grid(g,h,i,j,k,l, nrow =1 ,
                rel_widths = c(1.1,1,1,1,1,1))



p <- plot_grid(fig_broad, B1, fig_specific, B2, nrow = 4)


png("images/recount3-mouse-1.png",width = 1200, height = 1000)
print(p)
dev.off()



o <- fig2b("Drd2", "Favorite genes", "RNA-Seq  ")
r <- fig2b("Fos", " ", " ")
q <- fig2b("Grin1", " ", " ")
n <- fig2b("Mc4r", " ", " ")
s <- fig2b("Oxtr", " ", " ")
m <- fig2b("Pomc", " ", " ")

p2 <- plot_grid(o,r,q,n,s,m, nrow = 1)

png("images/recount3-mouse-2.png",width = 1200, height = 400)
print(p2)
dev.off()

# widen for tsne

fig_M1M2 <- png::readPNG("images/recount3-M1M2.png")
fig_M1M2 <- ggdraw() +  draw_image(fig_M1M2, scale = 1)



mycolors2 <- c( "#3b841e",  "#3d98b2",  "#0a1550" )

countData_long_wide <- countData_long %>%
  filter(!str_detect(mgi_symbol, "mt|Mir")) %>%
  filter(tissue %in% c("ca3d", "ca2", "ca1d"))  %>%
  dplyr::select(-mgi_symbol) %>%
  pivot_wider(id_cols = external_id:ensembl_gene_id,
              names_from = ensembl_gene_id, 
              values_from = counts,
              values_fn = sum)

tsne_data <- countData_long_wide[ ,10:50840]  %>%
  as.matrix(.)

tsne_samples <- countData_long_wide[ ,1:9] 
head(tsne_samples)


## Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_results <- Rtsne(tsne_data, perplexity=2, 
                      check_duplicates = FALSE) 

tsne_results_samples <- as.data.frame(tsne_results$Y) %>%
  cbind(tsne_samples, .)
head(tsne_results_samples) 

n <- tsne_results_samples %>%
  ggplot(aes(x = V1, y = V2, color = tissue)) +
  geom_point(size = 8, alpha = 0.75) +
  labs(x = "tSNE dimention 1", 
       y = "tSNE dimention 2", 
       color = "Replicates") +
  scale_color_manual(values = mycolors2,
                    labels = c("CA3", "CA2", "CA1")
                    ) +
  theme_linedraw(base_size = 14) 

p3 <- plot_grid(fig_M1M2, n, rel_widths = c(1, 1.25))
p3

png("images/recount3-mouse-3.png", width = 600, height = 250)
print(p3)
dev.off()
