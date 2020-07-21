#' R6 class for importing salmon count data
#'
#' @param project NgsProjectClass object
#' @param quantdir folder containing salmon quant data
#'
#' @return
#' @export
#'
#' @examples
#' project <- NgsProjectClass$new(projectdir)
SalmonClass <- R6::R6Class("SalmonClass",

  public = list(

    project = NULL,
    samples = NULL,
    quantdir = NULL,
    files = NULL,
    txi = NULL,

    initialize = function(project, quantdir)
    {
      print(paste0('quantdir=', quantdir))
      self$project <- project
      self$samples <- project$samples
      self$quantdir <- quantdir

      self$files <- private$listQuantFiles(project, quantdir)
      tx_exp <- tximport::tximport(self$files, type='salmon', txOut=TRUE) # dtuScaledTPM
      self$txi <- tximport::summarizeToGene(tx_exp, tx2gene=txdb$tx2gene, ignoreTxVersion=TRUE, countsFromAbundance='scaledTPM')

      invisible(self)
    },

    save = function(filename)
    {
      saveRDS(self, file = filename)
      invisible(self)
    },

    subset = function(samples=NULL, genes=NULL)
    {
      #print(paste0('samples=', samples, ' genes=', genes))
      #https://adv-r.hadley.nz/r6.html
      salmon <- self$clone(deep=TRUE)

      txi <- salmon$txi
      if (is.null(samples))
        samples <- salmon$samples$name
      if (is.null(genes))
        genes <- rownames(txi$counts)

      genes <- rownames(txi$counts)[rownames(txi$counts) %in% genes]

      txi$abundance <- txi$abundance[genes, samples]
      txi$counts <- txi$counts[genes, samples]
      txi$length <- txi$length[genes, samples]

      salmon$txi <- txi
      salmon$samples <- salmon$samples[samples,]
      return(salmon)
    }
  ),

  private = list(

    # make a list of the quant.sf files in quantdir
    listQuantFiles = function(project, quantdir)
    {
      # make a list of the quant.sf files under the current directory
      files <- file.path(quantdir, project$samples$name, "quant.sf")

      # confirm that the files exist
      all(file.exists(files))

      # tag each file with its sample name
      names(files) <- project$samples$name

      # check to make sure there is a corresponding quant.sf for each sample
      for (sample in project$samples$name)
      {
        if (!file.exists(files[[sample]]))
          R.oo::throw(concat(files[[sample]],' not found'))
      }

      return(files)
    }
  )
)
