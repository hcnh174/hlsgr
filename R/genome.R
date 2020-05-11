#data('genesymbol', package = "biovizBase")
#GENESYMBOL <- biovizBase::genesymbol
#hg19sub <- biovizBase::hg19sub


#' Split list of genes into a vector
#'
#' Splits comma-delimited list and returns only genes found in GENESYMBOL list
#'
#' @param values
#'
#' @return
#' @export
#'
#' @examples
#' hlsgr::splitGenes('ALK,NR1')#
#' hlsgr::splitGenes(c('PIK3CA', 'APC', 'FGFR2'))
splitGenes <- function(values=character(0))
{
  values <- splitFields(values)
  if (length(values)==0)
    return(character(0))
  return(values[values %in% names(biovizBase::genesymbol)])
}

#################################################################################

#' Create a circular chromsome plot showing location of variants
#'
#' https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
#'
#' @param snv
#' @param cnv
#' @param fusions
#'
#' @return
#' @export
#'
#' @examples
#' p <- hlsgr::plotChromosomes(cnv=c('GNA11', 'KRAS'), snv=c('PIK3CA', 'APC', 'FGFR2'), fusions=c('ALK', 'NBR1'))
plotChromosomes <- function(snv=character(0), cnv=character(0), fusions=character(0))
{
  snv <- splitGenes(snv)
  cnv <- splitGenes(cnv)
  fusions <- splitGenes(fusions)

  p <- ggbio::ggbio(trackWidth = 10, buffer = 0, radius = 10)
  if (length(fusions)>0)
    p <- p + ggbio::circle(biovizBase::genesymbol[fusions], geom = "rect", color = "green")
  if (length(cnv)>0)
    p <- p + ggbio::circle(biovizBase::genesymbol[cnv], geom = "rect", color = "red")
  if (length(snv)>0)
    p <- p + ggbio::circle(biovizBase::genesymbol[snv], geom = "rect", color = "orange")
  p <- p + ggbio::circle(biovizBase::hg19sub, geom = "ideo", fill = "gray70")
  p <- p + ggbio::circle(biovizBase::hg19sub, geom = "scale", size = 2)
  p <- p + ggbio::circle(biovizBase::hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
  return(p)
}


