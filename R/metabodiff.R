#' R6 class for MetaboDiff analysis
#'
#' @param samples samples
#'
#' @return
#' @export
#'
#' @examples
#' project <- MetaboDiffClass$new()
MetaboDiffClass <- R6::R6Class("MetaboDiffClass",

  inherit = AbstractBaseClass,

  public = list(

    samples =  NULL,

    initialize = function(samples)
    {
      super$initialize()
      self$samples <- samples
      invisible(self)
    }
  ),

  private = list(

  )
)
