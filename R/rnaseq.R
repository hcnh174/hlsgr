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

	return(list(salmon=salmon, dds=dds, res=res, res_lfc=res_lfc, mds=mds, vsd=vsd, annot_col=annot_col))
}

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
  genes <- rownames(result$res_lfc[order(abs(result$res_lfc$log2FoldChange), decreasing=TRUE), ])[1:num]
  pheatmap::pheatmap(SummarizedExperiment::assay(result$vsd)[genes,],
                     cluster_rows=TRUE, show_rownames=TRUE,
                     cluster_cols=TRUE, show_colnames=FALSE, scale='row',
                     annotation_col=result$annot_col, fontsize_row=fontsize,
                     color=color, clustering_distance_rows="correlation",
                     clustering_method='complete')
}

####################################################################################

#' Load sample information
#'
#' @param samples
#' @param countdir
#'
#' @return
#' @export
#'
#' @examples
#' filename <- 'samples.txt'
#' samples <- loadSamples(filename)
loadSamples <- function(filename)
{
  samples <- loadDataFrame(filename)
  samples$name <- as.character(samples$name)
  samples$subject <- as.character(samples$subject)
  samples$sample <- samples$name
  samples[samples==''] <- NA
  rownames(samples) <- samples$sample
  return(samples)
}

###################################################

#' loadSubjects
#'
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
#' groups <- loadSubjects('subjects.txt')
loadSubjects <- function(filename)
{
  subjects <- loadDataFrame(filename)
  return(subjects)
}

#########################################################

#' loadGroups
#'
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
#' groups <- loadGroups('treatments.txt')
loadGroups <- function(filename)
{
  groups <- loadDataFrame(filename)
  return(groups)
}

#####################################################

#' Load NGS project metadata from dir
#'
#' @param dir
#'
#' @return
#' @export
#'
#' @examples
loadNgsProject <- function(dir)
{
  subjects <- loadSubjects(paste0(dir, '/subjects.txt'))
  samples <- loadSamples(paste0(dir, '/samples.txt'))
  groups <- loadGroups(paste0(dir, '/groups.txt'))
  project <- list(samples = samples, groups=groups)
  return(project)
}

##############################################

#' Load Subread count files
#'
#' @param samples
#' @param countdir
#'
#' @return
#' @export
#'
#' @examples
#' countdir <-
#' subread <- loadSubreadData(samples, countdir)
loadSubreadData <- function(project, countdir)
{
  # make a list of the count files under the specified directory
  samples <- project$samples
  counts <- NA
  for (sample in samples$sample)
  {
    filename <- file.path(countdir, sample, paste0(sample, '.txt'))
    df <- as.data.frame(readr::read_tsv(filename, skip = 1)[, c(1,7)])
    colnames(df) <- c('geneid', as.character(sample))
    if (is.na(counts))
      counts <- df
    else counts <- cbind(counts, df[[sample]])
  }
  colnames(counts) <- c('geneid', samples$sample)
  rownames(counts) <- as.character(counts$geneid)
  counts <- counts[,-1]
  return(counts)
}

##################################################

#' getSamplesByGroup
#'
#' @param samples loaded from samples.txt
#' @param groups loaded from groups.txt
#' @param groupcol group column
#' @param excluded vector of samples to exclude
#'
#' @return
#' @export
#'
#' @examples
getSamplesByGroup <- function(project, groupcol, excluded=c())
{
  samples <- project$samples
  groups <- project$groups
  samples <- samples[!(samples$sample %in% exclude_samples),]
  samples$group <- samples[[paste0('group.', groupcol)]]
  samples <- samples[!is.na(samples$group),]
  samples$group <- factor(samples$group, levels=groups[groups$group==groupcol, 'level'])
  return(samples)
}
