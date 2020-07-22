#' create R6 class method for WGCNA
#'
#' @param result
#' @param datTraits
#'
#' @export
#'
#' @examples
#' wgcna <- WgcnaClass$new(result, datTraits)
WgcnaClass <- R6::R6Class("WgcnaClass",

  public = list(

    identifier = NULL,
    powers = c(c(1:10), seq(from = 12, to=20, by=2)),
    data = NULL,
    datTraits = NULL,
    gsg = NULL,
    matrix = NULL,
    sampleTree = NULL,
    cutHeight = NA,
    sft = NULL,
    net = NULL,
    moduleColors = NULL,
    hubgenes = NULL,

    initialize = function(rnaseq, datTraits, identifier='ensembl_id')
    {
      rnaseq <- private$checkMissing(rnaseq)

      wgcna_matrix <- t(rnaseq)
      print(dim(wgcna_matrix))
      sampleTree <- hclust(dist(wgcna_matrix), method = 'average')
      #datTraits <- private$prepareDatTraits(wgcna_matrix, project, group)

      self$identifier <- identifier
      self$data <- rnaseq
      self$matrix <- wgcna_matrix
      self$datTraits <- datTraits
      self$sampleTree <- sampleTree

      invisible(self)
    },

    plotSampleTree = function()
    {
      plot(self$sampleTree, main = "Sample clustering to detect outliers",
           sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
      #abline(h = 45, col = "red")
    },

    cutoffline = function(){
      cutHeight = readline("Enter cut line\n")
      cutHeight = as.numeric(unlist(strsplit(cutHeight, "[ ,]")))
      # Plot a line to show the cut
      abline(h = cutHeight, col = "red")
      self$cutHeight = cutHeight
      invisible(self)
    },

    cutTreeStatic = function()
    {
      clust <- WGCNA::cutreeStatic(self$sampleTree, cutHeight = self$cutHeight, minSize = 10)
      print(table(clust))
      keepSamples <- (clust==1)
      self$matrix <- self$matrix[keepSamples, ]
      self$datTraits <- self$datTraits[rownames(self$matrix),]
      invisible(self)
    },

    cutTreeDynamic = function(deepSplit=2)
    {
      clust <- dynamicTreeCut::cutreeDynamic(self$sampleTree, method='hybrid', deepSplit=deepSplit,
                                             distM=as.matrix(dist(self$matrix)), cutHeight=wgcna$sampleTree$height[length(wgcna$sampleTree$height)-2]-5)
      print(table(clust))
      keepSamples <- (clust!=0)
      self$matrix <- self$matrix[keepSamples, ]
      self$datTraits <- self$datTraits[rownames(self$matrix),]
      invisible(self)
    },

    plotDendogram = function(signed=FALSE)
    {
      ## Re-cluster samples
      sampleTree2 <- hclust(dist(self$matrix), method = "average")

      # Convert traits to a color representation: white means low, red means high, grey means missing entry
      traitColors = WGCNA::numbers2colors(self$datTraits, signed = signed)

      # Plot the sample dendrogram and the colors underneath.
      WGCNA::plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(self$datTraits),
                                 main = "Sample dendrogram and trait heatmap")
    },

    pickSoftThreshold = function()
    {
      sft <- WGCNA::pickSoftThreshold(self$matrix, powerVector = self$powers,
                                      RsquaredCut = 0.875, verbose = 5)
      self$sft <- sft
      invisible(self)
    },

    plotSoftThreshold = function()
    {
      sft <- self$sft
      plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           xlab='Soft Threshold (power)',
           ylab='Scale Free Topology Model Fit,signed R^2',
           type='n', main = paste('Scale independence'))
      text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           labels=self$powers, cex=1,col='red')
      abline(h=0.90,col='red')
    },

    plotMeanConnectivity = function()
    {
      sft <- self$sft
      plot(sft$fitIndices[,1], sft$fitIndices[,5],
           xlab="Soft Threshold (power)",
           ylab="Mean Connectivity", type="n",
           main = "Mean connectivity")
      text(sft$fitIndices[,1], sft$fitIndices[,5], labels=self$powers, cex=0.9,col="red")
    },

    constructNetwork = function(minModuleSize = 30, maxBlockSize = 5000, mergeCutHeight = 0.25, TOMType = "unsigned")
    {
      cor <- WGCNA::cor
      net <- WGCNA::blockwiseModules(self$matrix, weights = NULL,
                              power = self$sft$powerEstimate,
                              TOMType = TOMType, #TOMType = "signed",
                              minModuleSize = minModuleSize,
                              maxBlockSize = maxBlockSize,
                              reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = FALSE, #saveTOMFileBase = paste0(outdir,"/multiomics_ono_sensei_wgcna"),
                              verbose = 3, trapErrors=TRUE)
      #table(net$colors)
      self$net <- net
      self$moduleColors <- WGCNA::labels2colors(self$net$colors)
      invisible(self)
    },

    plotDendogramWithColors = function()
    {
      net <- self$net
      # Convert labels to colors for plotting
      #mergedColors <-  WGCNA::labels2colors(net$colors)

      # Plot the dendrogram and the module colors underneath
      WGCNA::plotDendroAndColors(net$dendrograms[[1]],
                                 self$moduleColors[net$blockGenes[[1]]],
                                 "Module colors", dendroLabels = FALSE,
                                 hang = 0.03, addGuide = TRUE, guideHang = 0.05)
    },

    plotTomPlot = function()
    {
      net <- self$net
      moduleLabels <- net$colors
      MEs <- net$MEs
      geneTree <- net$dendrograms[[1]]
      dissTOM <- 1 - WGCNA::TOMsimilarityFromExpr(self$matrix,
                                                  power = self$sft$powerEstimate, verbose = 3)
      plotTOM <- dissTOM^7
      diag(plotTOM) <- NA
      geneTree <- hclust(as.dist(dissTOM), method = "average")
      WGCNA::TOMplot(plotTOM, geneTree, self$moduleColors,
                     main = "Network heatmap plot, differentially expressed genes")
    },

    plotModuleTraitRelationship = function()
    {
      MEs0 <- moduleEigengenes(self$matrix, self$moduleColors)$eigengenes

      MEs <- orderMEs(MEs0)
      moduleTraitCor <- cor(MEs, self$datTraits, use = "p")
      moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(self$matrix))

      # Will display correlations and their p-values
      textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
      dim(textMatrix) <- dim(moduleTraitCor)

      if(dim(moduleTraitCor)[2] > 1)
      {
        # Display the correlation values within a heatmap plot
        WGCNA::labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(self$datTraits),
        yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE,
        colors = blueWhiteRed(50), colorMatrix = NULL, textMatrix = textMatrix,
        setStdMargins = TRUE, cex.text = 0.5, cex.lab = 0.7, zlim = c(-1,1),
        main = paste("Module-trait relationships"))
      }
      else
      {
        #only trait data is one column(eg. Control(0) and Case(1))
        image.plot(moduleTraitCor, xaxt="n", yaxt="n", col=colorRampPalette(c("blue", "white", "red"))(50),
        main = paste("Module-trait relationships"),
        breaks=seq(-max(moduleTraitCor), max(moduleTraitCor), 2*max(moduleTraitCor)/50))
        axis(side=1, at=seq(0,1, 1/(nrow(moduleTraitCor)-1)), labels=rownames(moduleTraitCor), las=2)
        axis(side=1, at=seq(0,1, 1/(nrow(moduleTraitCor)-1)), labels=textMatrix, font=2, pos=0.2, tick=FALSE)
        axis(side=2, at=0, labels="Control vs Case")
      }
    },

    selectHubs = function()
    {
      hub <- WGCNA::chooseTopHubInEachModule(self$matrix, self$moduleColors,
                                         power = self$sft$powerEstimate)
      hub <- data.frame(hub)
      colnames(hub) <- self$identifier
      annot <- txdb$gene[,c('ensembl_id', 'gene', 'description')]
      hubgenes <- merge(hub, annot, by=self$identifier, sort=FALSE)
      rownames(hubgenes) <- rownames(hub)
      self$hubgenes <- hubgenes
      outfile <- paste0(outdir, '/hubgenes.txt')
      writeTable(self$hubgenes, outfile, row.names=TRUE)
      invisible(self)
    },

    goEnrichmentAnalysis = function()
    {
      annot <- txdb$gene[,c("ensembl_id", "entrez_id", "gene")]
      # Match probes in the data set to the probe IDs in the annotation file
      probes <- colnames(self$matrix)
      probes2annot <-  match(probes, annot[[self$identifier]])
      # Get the corresponding Locuis Link IDs
      allLLIDs <- annot$entrez_id[probes2annot]
      # Choose interesting modules
      #intModules = c("brown", "red", "salmon")
      intModules <- unique(self$moduleColors[self$moduleColors != "grey"])

      ###Enrichment analysis directly within R
      go <- WGCNA::GOenrichmentAnalysis(self$moduleColors, allLLIDs, organism = "human", nBestP = 10)
      tab <- go$bestPTerms[[4]]$enrichment
      screenTab <- tab[, c('module', 'modSize', 'enrichmentP', 'BonferoniP', 'nModGenesInTerm', 'termOntology', 'termName')]
      # Round the numeric columns to 2 decimal places:
      numericColumns <- c('enrichmentP', 'BonferoniP')
      screenTab[, numericColumns] <- signif(apply(screenTab[, numericColumns], 2, as.numeric), 2)
      # Truncate the the term name to at most 40 characters
      screenTab[, 'termName'] <- substring(screenTab[, 'termName'], 1, 40)
      # Shorten the column names:
      colnames(screenTab) <- c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name")
      rownames(screenTab) <- NULL
      outfile <- paste0(outdir, '/go_enrichment.txt')
      writeTable(screenTab, outfile)
      return(screenTab)
    },

    #Export Cytoscape
    exportCytoscape = function(threshold = 0.02)
    {
      # Recalculate topological overlap if needed
      TOM <- WGCNA::TOMsimilarityFromExpr(self$matrix, power = self$sft$powerEstimate, verbose = 3)
      # Read in the annotation file
      annot <- txdb$gene[,c("ensembl_id","gene")]
      colors <- unique(self$moduleColors)
      unique_colors <- colors[colors != "grey"]
      for (module in unique_colors)
      {
        # Select module probes
        probes <- colnames(self$matrix)
        inModule <- is.finite(match(self$moduleColors, module))
        modProbes <- probes[inModule]
        modGenes <- annot$gene[match(modProbes, annot[[self$identifier]])]

        # Select the corresponding Topological Overlap
        modTOM <- TOM[inModule, inModule]
        dimnames(modTOM) <- list(modProbes, modProbes)

        # Export the network into edge and node list files Cytoscape can read
        edgeFile <- paste0(outdir, '/CytoscapeInput-edges-', module, '.txt')
        nodeFile <- paste0(outdir, '/CytoscapeInput-nodes-', module, '.txt')
        cyt <- exportNetworkToCytoscape(modTOM, edgeFile = edgeFile, nodeFile = nodeFile,
          weighted = TRUE, threshold = threshold, nodeNames = modProbes,
          altNodeNames = modGenes, nodeAttr = self$moduleColors[inModule])
      }
    }
  ),

  private = list(

    checkMissing = function(datExpr)
    {
      gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
      print(gsg$allOK)
      if (!gsg$allOK)
      {
        # Optionally, print the gene and sample names that were removed:
        if (sum(!gsg$goodGenes)>0)
          printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
        if (sum(!gsg$goodSamples)>0)
          printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
        # Remove the offending genes and samples from the data:
        datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
      }
      self$gsg <- gsg
      return(datExpr)
    }

    # prepareDatTraits = function(wgcna_matrix, project, group)
    # {
    #   samples <- hlsgr::getSamplesByGroup(project, group)
    #   samples <- samples[rownames(samples) %in% rownames(wgcna_matrix), ]
    #   return(samples)
    # }
  )
)

