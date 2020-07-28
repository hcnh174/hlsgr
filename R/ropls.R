#' R6 class for RoplsClass analysis
#'
#' @param samples samples
#'
#' @return
#' @export
#'
#' @examples
#' project <- RoplsClass$new()
RoplsClass <- R6::R6Class("RoplsClass",

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
