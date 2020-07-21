#library(R6)
#library(mixOmics)

#' R6 class for mixomics analysis
#'
#' @param result DGE analysis result
#'
#' @return
#' @export
#'
#' @examples
#' mixomixs <- MixomicsClass$new(result)
MixomicsClass <- R6::R6Class("MixomicsClass",
  public = list(

    data = NULL,
    samples = NULL,
    design = NULL,
    ncomp = NULL,
    sgccda.res = NULL,
    perf.diablo = NULL,
    list.keepX = NULL,

    initialize = function(result) {
      count_data <- DESeq2::counts(result$dds, normalized=TRUE)
      keep <- rowSums(counts(result$dds)) >= 70
      count_data <- count_data[keep,]

      samples <- result$salmon$samples
      samples <- samples[!is.na(samples$group),]

      max_genes <- 10000
      mRNA <- t(count_data[order(apply(count_data, 1, mad), decreasing = T)[1:max_genes],])

      metabolomics <- read.csv(paste0(metabolomicsdir, "/metabolome_normalization.txt"), header=T, check.names=F, row.names=1)

      multi <- list(
        mRNA = mRNA,
        metabolomics = t(metabolomics))

      lapply(multi, dim)

      ids <- base::Reduce(intersect, list(rownames(multi$mRNA),
                                          rownames(multi$metabolomics)))

      multi$mRNA <- multi$mRNA[ids, ]
      multi$metabolomics <- multi$metabolomics[ids, ]

      lapply(multi, dim)

      samples <- samples[ids,]

      samples$group <- factor(samples$group, levels=c("Control", "Case"))

      self$data <- multi
      self$samples <- samples
      self$design <- private$createDesign(multi)

      invisible(self)
    },

    evaluatePerformance = function()
    {
      self$sgccda.res <- mixOmics::block.splsda(X = self$data,
                                                Y = self$samples$group,
                                                ncomp = 5, design = self$design)

      # this code takes a few minutes to run
      self$perf.diablo <- mixOmics::perf(mixdata$sgccda.res, validation = 'Mfold',
                                         folds = 10, nrepeat = 10)
    },

    plotPerformance = function()
    {
      plot(self$perf.diablo)
    },

    tune = function()
    {
      # $choice.ncomp indicates an optimal number of components for the final DIABLO model.
      #self$perf.diablo$choice.ncomp$WeightedVote

      # optimal number of components
      self$ncomp <- self$perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

      # model tuning
      #test.keepX = list(mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
      #		metabolomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))

      test.keepX = list (mRNA = c(5, 10, 15, 20),
                         metabolomics = c(5, 10, 15, 20))

      tune.TCGA = mixOmics::tune.block.splsda(X = self$data, Y = self$samples$group,
              ncomp = self$ncomp, test.keepX = test.keepX, design = self$design,
              validation = 'Mfold', folds = 10, nrepeat = 1, cpus = 2,
              dist = "centroids.dist")

      self$list.keepX <- tune.TCGA$choice.keepX
      #list.keepX = list(mRNA = c(9, 6, 6), metabolomics = c(30, 7, 5))
      print(self$list.keepX)
      #mRNA: 15  5 20 20, metabolomics: 10 20 10
      invisible(self)
    },

    getFinalModel = function()
    {
      self$sgccda.res <- mixOmics::block.splsda(X = self$data,
        Y = self$samples$group, ncomp = self$ncomp,
        keepX = self$list.keepX, design = self$design)
      #self$sgccda.res$design
      invisible(self)
    },

    plotDiablo = function(ncomp = 1)
    {
      mixOmics::plotDiablo(self$sgccda.res, ncomp = ncomp)
      invisible(self)
    },

    plotIndiv = function()
    {
      mixOmics::plotIndiv(self$sgccda.res, ind.names = FALSE,
                          legend = TRUE, title = 'DIABLO')
      invisible(self)
    },

    plotArrow = function()
    {
      mixOmics::plotArrow(self$sgccda.res, ind.names = FALSE,
                          legend = TRUE, title = 'DIABLO')
      invisible(self)
    },

    plotVar = function()
    {
      mixOmics::plotVar(self$sgccda.res, var.names = FALSE, style = 'graphics',
                        legend = TRUE, pch = c(16, 17), cex = c(2,2),
                        col = c('darkorchid', 'lightgreen'))
      invisible(self)
    },

    plotCircos = function()
    {
      sgccda.res <- private$annotateResults(self$sgccda.res)
      mixOmics::circosPlot(sgccda.res, cutoff = 0.4, line = TRUE,
                           color.blocks= c('darkorchid','lightgreen'),
                           color.cor = c("chocolate3","grey20"),
                           size.labels = 1.5, size.variables = 0.5)
      invisible(self)
    },

    plotLoadings = function()
    {
      sgccda.res <- private$annotateResults(self$sgccda.res)
      for(i in 1:self$ncomp)
      {
        mixOmics::plotLoadings(sgccda.res, comp = i, contrib = 'max',
          method = 'median', size.name = 1.0)#, col=self$col
      }
      invisible(self)
    },

    plotNetwork = function()
    {
      sgccda.res <- private$annotateResults(self$sgccda.res)
      net <- mixOmics::network(sgccda.res, blocks = c(1,2),
                               color.node = c('darkorchid', 'lightgreen'), cutoff = 0.4)
      #write.graph(net$gR, file = paste(outdir, "myNetwork.gml", sep=""), format = "gml")
      invisible(self)
    },

    plotHeatmap = function()
    {
      sgccda.res <- private$annotateResults(self$sgccda.res)
      mixOmics::cimDiablo(sgccda.res, margins=c(5,20), color.blocks=c('pink', 'green'))
      invisible(self)
    }
  ),

  private = list(
    createDesign = function(data)
    {
      design <- matrix(0.1, ncol = length(data), nrow = length(data),
                       dimnames = list(names(data), names(data)))
      diag(design) = 0
      return(design)
    },

    annotateResults = function(sgccda.res.orig)
    {
      sgccda.res <- sgccda.res.orig

      ensembl_ids <- rownames(sgccda.res$loadings$mRNA)
      genes <-  lookupEnsemblIds(ensembl_ids)
      rownames(sgccda.res$loadings$mRNA) <- genes
      colnames(sgccda.res$X$mRNA) <- genes
      sgccda.res$names$colnames$mRNA <- genes

      hmt_ids <- rownames(sgccda.res$loading$metabolomics)
      hmt_names <- lookupHmtIds(hmt_ids)
      rownames(sgccda.res$loadings$metabolomics) <- hmt_names
      colnames(sgccda.res$X$metabolomics) <- hmt_names
      sgccda.res$names$colnames$metabolomics <- hmt_names

      return(sgccda.res)
    }
  )
)

