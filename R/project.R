#' R6 class for project information
#'
#' @param dir folder containing project and sample information
#'
#' @return
#' @export
#'
#' @examples
#' project <- NgsProjectClass$new(projectdir)
NgsProjectClass <- R6::R6Class("NgsProjectClass",

   public = list(

     subjects = NULL,
     samples = NULL,
     runs = NULL,
     groups = NULL,

     initialize = function(dir)
     {
       print(paste0('dir=', dir))
       self$subjects <- private$loadSubjects(paste0(dir, '/subjects.txt'))
       self$samples <- private$loadSamples(paste0(dir, '/samples.txt'))
       self$runs <- private$loadRuns(paste0(dir, '/runs.txt'))
       self$groups <- private$loadGroups(paste0(dir, '/groups.txt'))
       invisible(self)
     },

     #' getGroupLevels returns unique ordered levels for specified group
     #'
     #' @param groupcol group column
     #'
     #' @return
     #' @export
     #'
     #' @examples
     #' getGroupLevels(project, groupcol='default')
     getGroupLevels = function(groupcol)
     {
        rows <- self$groups[self$groups$group==groupcol,]
        if (nrow(rows)==0)
           R.oo::throw('groupcol: ', groupcol, ' not found in groups dataframe')
        return(rows[['level']])
     },

     ##################################################

    #' getSamplesByGroup
    #'
    #' @param groupcol group column
    #'
    #' @return
    #' @export
    #'
    #' @examples
    getSamplesByGroup = function(groupcol)
    {
      samples <- self$samples
      #groups <- project$groups
      #samples <- samples[!(samples$sample %in% excluded),]
      samples$group <- samples[[paste0('group.', groupcol)]]
      samples <- samples[!is.na(samples$group),]
      samples$group <- factor(samples$group, levels=self$getGroupLevels(groupcol))
      samples$groupnum <- as.numeric(samples$group) -1
      return(samples)
    },

    getCountsByGroup = function(counts, groupcol)
    {
      # subset samples non-null for the selected groupcol
      samples <- self$getSamplesByGroup(groupcol)
      # susbet the count data to include only samples selected above
      counts <- counts[, colnames(counts) %in% samples$sample]
      # makes sure that the sample list also matches samples in the count data
      samples <- samples[samples$name %in% colnames(counts),]

      counts <- counts[, rownames(samples)]
      print(all(rownames(samples) %in% colnames(counts)))
      print(all(rownames(samples) == colnames(counts)))

      return(list(samples = samples, counts = counts, groupcol = groupcol))
    }
   ),

   private = list(

     loadSamples = function(filename)
     {
       samples <- loadDataFrame(filename)
       samples$name <- as.character(samples$name)
       samples$subject <- as.character(samples$subject)
       samples$sample <- samples$name
       samples[samples==''] <- NA
       rownames(samples) <- samples$sample
       return(samples)
     },

     ###################################################

     loadSubjects = function(filename)
     {
       subjects <- loadDataFrame(filename)
       return(subjects)
     },

     #########################################################

     loadGroups = function(filename)
     {
       groups <- loadDataFrame(filename)
       return(groups)
     },

     #########################################################

     loadRuns = function(filename)
     {
       runs <- loadDataFrame(filename)
       return(runs)
     }
   )
)

#######################################################

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
#' subjects <- loadSubjects('subjects.txt')
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

#########################################################

#' loadRuns
#'
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
#' groups <- loadRuns('runs.txt')
loadRuns <- function(filename)
{
  runs <- loadDataFrame(filename)
  return(runs)
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
  runs <- loadRuns(paste0(dir, '/runs.txt'))
  groups <- loadGroups(paste0(dir, '/groups.txt'))
  project <- list(subjects = subjects, samples = samples, runs = runs, groups=groups)
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

#' getGroupLevels returns uniquq ordered levels for specified group
#'
#' @param project project object
#' @param groupcol group column
#'
#' @return
#' @export
#'
#' @examples
#' getGroupLevels(project, groupcol='default')
getGroupLevels <- function(project, groupcol)
{
  return(project$groups[project$groups$group==groupcol, 'level'])
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
  #groups <- project$groups
  samples <- samples[!(samples$sample %in% excluded),]
  samples$group <- samples[[paste0('group.', groupcol)]]
  samples <- samples[!is.na(samples$group),]
  samples$group <- factor(samples$group, levels=getGroupLevels(project, groupcol))
  samples$groupnum <- as.numeric(samples$group) -1
  return(samples)
}

##################################################################

#' getCountsByGroup adjusts the counts and samples to exclude missing samples
#'
#' @param counts RNAseq count data
#' @param project project data containing sample and group information
#' @param groupcol group column
#'
#' @return
#' @export
#'
#' @examples
getCountsByGroup <- function(counts, project, groupcol)
{
  # subset samples non-null for the selected groupcol
  samples <- hlsgr::getSamplesByGroup(project, groupcol)
  # susbet the count data to include only samples selected above
  counts <- counts[, colnames(counts) %in% samples$sample]
  # makes sure that the sample list also matches samples in the count data
  samples <- samples[samples$name %in% colnames(counts),]
  return(list(samples = samples, counts = counts, groupcol = groupcol))
}
