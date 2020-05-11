# load protein-coding genes
# Question: List of Ensembl Transcript IDs corresponding to protein-coding genes
# https://www.biostars.org/p/173427/
protein_coding_genes <- hlsgr::loadDataFrame(paste0('data-raw/protein_coding_genes.txt'))
#include_genes <- coding$gene_id

usethis::use_data(protein_coding_genes, overwrite = TRUE)
