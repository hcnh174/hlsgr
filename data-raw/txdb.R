## code to prepare `DATASET` dataset goes here

#' Load data table linking Ensmebl transcript and gene IDs
#'
#' @param quantdir
#'
#' @return
#' @export
#'
#' @examples
#' txdb <- loadTxdb()
loadTxdb <- function()
{
  txdb <- loadDataFrame(concat('data-raw/tx2gene.txt'), quote='')
  tx2gene <- txdb[,c('ensembl_transcript_id', 'ensembl_gene_id')]
  colnames(tx2gene) <- c('TXNAME', 'GENEID')
  #print(head(tx2gene))

  # also map ensemble_gene_id to entrez_id
  genes <- txdb[,c('ensembl_gene_id', 'entrezgene', 'external_gene_name', 'description')]
  genes <- unique(genes)
  colnames(genes) <- c('nsembl_id', 'entrez_id', 'gene', 'description')
  genes$description <- gsub(' \\[..+', '', genes$description)
  #print(head(genes))

  return(list(txdb=txdb, tx2gene=tx2gene, genes=genes))
}

txdb <- loadTxdb()

usethis::use_data(txdb, overwrite = TRUE)

