# library(ggbio)
# library(edgeR)
# library(limma)
# library(csaw)

#https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
#https://rdrr.io/bioc/edgeR/man/catchSalmon.html
#https://angus.readthedocs.io/en/2017/counting.html

#' R6 class for EdgeR analysis
#'
#' @param salmon salmon count data
#' @param groupcol category to use for DGE contrasts
#' @param outdir directory to write DGE results to
#'
#' @return
#' @export
#'
#' @examples
#' project <- EdgeRClass$new(salmon, groupcol, outdir)
EdgeRClass <- R6::R6Class("EdgeRClass",

  public = list(

    salmon = NULL,
    samples = NULL,
    design = NULL,
    fit = NULL,

    initialize = function(salmon, groupcol, outdir)
    {
      #https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#edgeR
      samples <- salmon$project$getSamplesByGroup(groupcol)
      salmon <- salmon$subset(samples=samples$name)

      #files <- stringr::str_replace(salmon$files, '/quant.sf', '')
      #quants <- edgeR::catchSalmon(files)
      #dge <- DGEList(counts=quants$counts/quants$annotation$Overdispersion, genes=quants$annotation)

      # cts <- salmon$txi$counts
      # normMat <- salmon$txi$length
      #
      # # Obtaining per-observation scaling factors for length, adjusted to avoid
      # # changing the magnitude of the counts.
      # normMat <- normMat/exp(rowMeans(log(normMat)))
      # normCts <- cts/normMat
      #
      # # Computing effective library sizes from scaled counts, to account for
      # # composition biases between samples.
      # eff.lib <- edgeR::calcNormFactors(normCts) * colSums(normCts)
      #
      # # Combining effective library sizes with the length factors, and calculating
      # # offsets for a log-link GLM.
      # normMat <- sweep(normMat, 2, eff.lib, "*")
      # normMat <- log(normMat)
      #
      # # Creating a DGEList object for use in edgeR.
      # y <- DGEList(cts)
      # y <- scaleOffset(y, normMat)
      # # filtering
      # keep <- filterByExpr(y)
      # ## Warning in filterByExpr.DGEList(y): All samples appear to belong to the same
      # ## group.
      # y <- y[keep, ]
      #
      # # y is now ready for estimate dispersion functions see edgeR User's Guide
      # #For creating a matrix of CPMs within edgeR, the following code chunk can be used:
      # se <- SummarizedExperiment(assays = list(counts = y$counts, offset = y$offset))
      # se$totals <- y$samples$lib.size
      #
      # cpms <- csaw::calculateCPM(se, use.offsets = TRUE, log = FALSE)
      #dge <- edgeR::DGEList(counts=cpms, group=samples$group)

      # https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
      dge <- edgeR::DGEList(counts=salmon$txi$counts, group=samples$group)

      # create a design matrix
      #design <- model.matrix(~ 0+samples$group)
      #colnames(design) <- levels(samples$group)
      design <- model.matrix(~ samples$group)
      colnames(design) <- levels(samples$group)

      # The next step is to remove rows that consistently have zero or very low counts
      keep <- edgeR::filterByExpr(dge, design)
      dge <- dge[keep, , keep.lib.sizes=FALSE]

      # limma-trend
      logCPM <- cpm(dge, log=TRUE, prior.count=3)
      #boxplot(logCPM)

      fit <- limma::lmFit(logCPM, design)
      fit <- limma::eBayes(fit, trend=TRUE)

      tbl <- limma::topTable(fit, coef=ncol(design), sort.by='logFC', number=500)
      groupname <- tolower(substring(colnames(design)[2], 14))
      outfile <- paste0(outdir, '/table-edger-group-',groupname,'.txt')
      hlsgr::writeTable(tbl, outfile, row.names=TRUE)

      self$salmon <- salmon
      self$samples <- samples
      self$design <- design
      self$fit <- fit
      invisible(self)
    },

    getResults = function(p=NULL, padj=0.05, logfc=1.5)
    {
      results <- limma::topTable(self$fit, coef=ncol(self$design), sort.by='logFC')
      results <- na.omit(results)
      if (!is.null(p))
        results <- results[results$P.Value < p, ]
      if (!is.null(padj))
        results <- results[results$adj.P.Val < padj, ]
      if (!is.null(logfc))
        results <- results[abs(results$logFC) > logfc, ]
      return(results)
    }
  ),

  private = list(

  )
)
