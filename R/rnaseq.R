#library(openxlsx)
#library(tximport)
#library(edgeR)
#library(DESeq2)
#library(pheatmap)

#' Load sample information for Salmon-mapped samples
#'
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
#' samples <- loadSalmonSamples()
loadSalmonSamples <- function(filename)#='data/samples.txt'
{
	samples <- loadDataFrame(filename)
	samples$sample <- as.character(samples$sample)
	samples$subject <- as.character(samples$subject)
	rownames(samples) <- samples$sample
	return(samples)
}

####################################################################################

#' Load Salmon quant files
#'
#' @param samples
#' @param quantdir
#'
#' @return
#' @export
#'
#' @examples
#' salmon <- loadSalmonData(samples, quantdir)
loadSalmonData <- function(samples, quantdir)
{
	#txdb <- loadTxdb(quantdir)

	# make a list of the quant.sf files under the current directory
	files <- file.path(quantdir, samples$folder, "quant.sf")

	# tag each file with its sample name
	names(files) <- samples$sample

	# check to make sure there is a corresponding quant.sf for each sample
	for (sample in samples$sample)
	{
		if (!file.exists(files[[sample]]))
			R.oo::throw(concat(files[[sample]],' not found'))
	}

	# import the quant data from salmon using tximport
	tx.exp <- tximport::tximport(files, type='salmon', txOut=TRUE) # dtuScaledTPM
	gene.exp <- tximport::summarizeToGene(tx.exp, tx2gene=txdb$tx2gene, ignoreTxVersion=TRUE, countsFromAbundance='scaledTPM')

	return(list(txi=gene.exp, samples=samples, txdb=txdb))
}

##############################################################################

#' Write Salmon quant data to to rds file for faster loading across sessions
#'
#' @param salmon salmon quant data
#' @param filename .rds file to store quant data in
#'
#' @return
#' @export
#'
#' @examples
#' saveSalmonData(salmon, 'salmon.rds')
saveSalmonData <- function(salmon, filename)
{
  saveRDS(salmon, file = filename)
}

############################################################

#' Read Salmon data from RDS file
#'
#' @param filename .rds file holding salmon quant data
#'
#' @return
#' @export
#'
#' @examples
#' salmon = readSalmonData('salmon.rds')
readSalmonData <- function(filename)
{
  return(readRDS(file = filename))
}

####################################################################################

#' Subset txi data by sample or gene
#'
#' @param txi
#' @param samples
#' @param include_genes
#'
#' @return
#' @export
#'
#' @examples
subsetTxi <- function(txi, samples, include_genes=rownames(txi$counts))
{
	genes <- rownames(txi$counts)[rownames(txi$counts) %in% include_genes]
	txi$abundance <- txi$abundance[genes, samples$sample]
	txi$counts <- txi$counts[genes, samples$sample]
	txi$length <- txi$length[genes, samples$sample]
	return(txi)
}

####################################################################################

#' Analyze Salmon quant data using DeSeq2
#'
#' @param salmon
#' @param exclude_samples
#' @param include_genes
#' @param groupcol
#' @param levels
#'
#' @return
#' @export
#'
#' @examples
analyzeSalmonDataDeSeq2 <- function(salmon, outdir, exclude_samples=c(), protein_only = FALSE,
                                    groupcol='group', levels=c('Control', 'Case'))
{
  print(paste0('outdir=', outdir))
  print(paste0('exclude_samples=', exclude_samples))
  print(paste0('protein_only=', protein_only))
  print(paste0('groupcol=', groupcol))
  print(paste0('levels=', levels))

	samples <- salmon$samples
	samples <- samples[!(samples$sample %in% exclude_samples),]
	samples <- samples[!is.na(samples[[groupcol]]),]
	samples[[groupcol]] <- factor(samples[[groupcol]], levels=levels)

	include_genes <- if(protein_only) protein_coding_genes$gene_id else rownames(salmon$txi$counts)

	txi <- subsetTxi(salmon$txi, samples, include_genes=include_genes)

	# create the design matrix
	group <- factor(samples[colnames(txi$counts), groupcol], levels=levels)
	design <- model.matrix(~group)
	print(design)

	# https://angus.readthedocs.io/en/2019/diff-ex-and-viz.html
	dds <- DESeq2::DESeqDataSetFromTximport(txi = txi, colData = samples, design = design)
	dds <- DESeq2::DESeq(dds)
	res <- DESeq2::results(dds)
	res <- res[order(res$log2FoldChange, decreasing=TRUE), ]

	groupname <- tolower(substring(colnames(design)[2], 6))
	outfile <- concat(outdir, '/table-deseq2-group-',groupname,'.txt')
	writeTable(as.data.frame(res), outfile, row.names=TRUE)

	res_sig <- subset(res, padj<.05)
	res_lfc <- subset(res_sig, abs(log2FoldChange) > 1)
	#res_sorted_abs <- res_lfc[order(abs(res_lfc$log2FoldChange), decreasing=TRUE), ]
	#head(res_lfc)

	vsd <- DESeq2::vst(dds)

	sample_dists <- SummarizedExperiment::assay(vsd) %>%
	  t() %>%
	  dist() %>%
	  as.matrix()

	mdsData <- data.frame(cmdscale(sample_dists))
	mds <- cbind(mdsData, as.data.frame(SummarizedExperiment::colData(vsd))) # combine with sample data

	annot_col <- samples[,c('sample', 'group')] %>%
	  dplyr::select(group) %>%
	  as.data.frame()

	return(list(salmon=salmon, samples=samples, dds=dds, res=res, res_lfc=res_lfc, mds=mds, vsd=vsd, annot_col=annot_col))
}

################################################

#' Analyze Salmon quant data using DeSeq2
#'
#' @param project
#' @param txi
#' @param groupcol
#' @param outdir
#' @param protein_only
#'
#' @return
#' @export
#'
#' @examples
analyzeSalmonDeSeq2 <- function(project, txi, groupcol, outdir, protein_only = FALSE)
{
  print(paste0('groupcol=', groupcol))
  print(paste0('outdir=', outdir))
  print(paste0('protein_only=', protein_only))

  samples <- getSamplesByGroup(project, groupcol)
  levels <- getGroupLevels(project, groupcol)
  print(paste0('levels=', levels))

  #samples <- samples[!(samples$sample %in% exclude_samples),]
  #samples <- samples[!is.na(samples[[groupcol]]),]
  #samples[[groupcol]] <- factor(samples[[groupcol]], levels=levels)

  include_genes <- if(protein_only) protein_coding_genes$gene_id else rownames(txi$counts)

  txi <- subsetTxi(txi, samples, include_genes=include_genes)

  # create the design matrix
  group <- factor(samples[colnames(txi$counts), 'group'], levels=levels)
  design <- model.matrix(~group)
  print(design)

  # https://angus.readthedocs.io/en/2019/diff-ex-and-viz.html
  dds <- DESeq2::DESeqDataSetFromTximport(txi = txi, colData = samples, design = design)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  res <- res[order(res$log2FoldChange, decreasing=TRUE), ]

  groupname <- tolower(substring(colnames(design)[2], 6))
  outfile <- concat(outdir, '/table-deseq2-group-',groupname,'.txt')
  writeTable(as.data.frame(res), outfile, row.names=TRUE)

  res_sig <- subset(res, padj<.05)
  res_lfc <- subset(res_sig, abs(log2FoldChange) > 1)
  #res_sorted_abs <- res_lfc[order(abs(res_lfc$log2FoldChange), decreasing=TRUE), ]
  #head(res_lfc)

  vsd <- DESeq2::vst(dds)

  sample_dists <- SummarizedExperiment::assay(vsd) %>%
    t() %>%
    dist() %>%
    as.matrix()

  mdsData <- data.frame(cmdscale(sample_dists))
  mds <- cbind(mdsData, as.data.frame(SummarizedExperiment::colData(vsd))) # combine with sample data

  annot_col <- samples[,c('sample', 'group')] %>%
    dplyr::select(group) %>%
    as.data.frame()

  return(list(txi=txi, samples=samples, dds=dds, res=res, res_lfc=res_lfc, mds=mds, vsd=vsd, annot_col=annot_col))
}

#####################################################

#' Create a heatmap for DeSeq2 RNASeq data
#'
#' @param result
#' @param num
#' @param fontsize
#' @param color
#'
#' @return
#' @export
#'
#' @examples
#' plotHeatmap(result, 10)
plotDeSeq2Heatmap <- function(result, num=50, fontsize=6, color=colorRampPalette(c("blue", "white", "red"))(50))
{
  ensembl_ids <- rownames(result$res_lfc[order(abs(result$res_lfc$log2FoldChange), decreasing=TRUE), ])[1:num]
  data <- SummarizedExperiment::assay(result$vsd)[ensembl_ids,]
  rownames(data) <- lookupGenes(ensembl_ids)
  pheatmap::pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE,
                     cluster_cols=TRUE, show_colnames=FALSE, scale='row',
                     annotation_col=result$annot_col, fontsize_row=fontsize,
                     color=color, clustering_distance_rows="correlation",
                     clustering_method='complete')
}

########################################

#' Return a vector gene names corresponding to a vector of Ensembl IDs
#'
#' @param ensembl_ids vector of ensembl IDs
#'
#' @return
#' @export
#'
#' @examples
#' lookupGenes(c('ENSG00000211459', 'ENSG00000211958 '))
lookupGenes <- function(ensembl_ids)
{
  genedf <- unique(txdb$genes[which(txdb$genes$ensembl_id %in% ensembl_ids), c('ensembl_id', 'gene')])
  rownames(genedf) <- genedf$ensembl_id
  return(genedf[ensembl_ids, 'gene'])
}

####################################################################################

#' Create a formatted table of differentially regulated genes
#'
#' @param result DeSeq2 analysis results
#' @param num max results to show
#' @param padj cutoff value for adjusted p-value
#' @param lfc log fold change cutoff (positive or negative)
#'
#' @return
#' @export
#'
#' @examples
#' makeDeSeq2Table(result)
makeDeSeq2Table <- function(result, num=50, padj=1, lfc=1.5)
{
  tbl <- result$res_lfc
  tbl <- tbl[tbl$padj<padj, ]

  if (lfc > 0)
  {
    tbl <- tbl[tbl$log2FoldChange>lfc, ]
    tbl <- tbl[order(tbl$log2FoldChange, decreasing=TRUE), ]
  }
  if (lfc < 0)
  {
    tbl <- tbl[tbl$log2FoldChange<lfc, ]
    tbl <- tbl[order(tbl$log2FoldChange), ]
  }

  tbl$ensembl_id <- rownames(tbl)

  tbl <- merge(as.data.frame(tbl), txdb$genes, by.x='ensembl_id', by.y='ensembl_id')
  tbl <- tbl[, c("gene", 'description', "log2FoldChange", "pvalue", "padj" )]

  tbl <- tbl[!duplicated(tbl), ]

  tbl <- head(tbl, n=num)
  return(tbl)
}

