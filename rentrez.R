#library(devtools)
#install_github("ropensci/rentrez")

library(rentrez)

paper <- entrez_search(db="pubmed", 
                       term="10.1016/j.ygcen.2013.10.007[doi]")
paper$ids

data <- entrez_link(db="all", id=paper$ids, dbfrom="pubmed")
data$links

proteins <- entrez_fetch(db="protein", 
                         id=data$links$pubmed_protein, 
                         rettype="fasta", )
cat(proteins)
