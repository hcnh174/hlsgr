#' R6 class for MetaboAnalyst analysis
#'
#' @param samples samples
#'
#' @return
#' @export
#'
#' @examples
#' project <- MetaboAnalystClass$new()
MetaboAnalystClass <- R6::R6Class("MetaboAnalystClass",

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
