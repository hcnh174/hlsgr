#' R6 class to provide a common interface for count data
#'
#' @param dir folder containing project and sample information
#'
#' @return
#' @export
#'
#' @examples
#' project <- RnaSeqCountDataClass$new(projectdir)
RnaSeqCountDataClass <- R6::R6Class("RnaSeqCountDataClass",

  inherit = AbstractBaseClass,

  public = list(

    counts = NULL,
    identifier = NULL,

    initialize = function(counts, identifier = 'ensembl_id')
    {
      self$counts = counts
      self$identifier = identifier
      invisible(self)
    }
  )
)
