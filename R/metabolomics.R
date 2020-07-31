
AbstractMetabolomicsDataClass <- R6::R6Class("AbstractMetabolomicsDataClass",

  inherit = AbstractBaseClass,

  public = list(
    filename = NULL,
    data = NULL,

    initialize = function(filename)
    {
      super$initialize()
      print('AbstractMetabolomicsDataClass')
      self$filename <- filename
      invisible(self)
    }
  )
)

#############################################

#' R6 class for Takeda metabolomics data
#'
#' @param filename filename
#'
#' @return
#' @export
#'
#' @examples
#' data_metabolomics <- TakedaMetabolomicsDataClass$new(paste0(dir, '/data.txt'))
TakedaMetabolomicsDataClass <- R6::R6Class("TakedaMetabolomicsDataClass",

  inherit = AbstractMetabolomicsDataClass,

  public = list(

    initialize = function(filename)
    {
      super$initialize(filename)
      print('TakedaMetabolomicsDataClass')
      self$data <- hlsgr::loadDataFrame(filename, check.names = FALSE)
	    rownames(self$data) <- self$data$id
      invisible(self)
    }
  )
)

#######################################################

#' R6 class for HMT metabolomics data
#'
#' @param filename filename
#'
#' @return
#' @export
#'
#' @examples
#' data_metabolomics <- HmtMetabolomicsDataClass$new(paste0(metabolomicsdir, '/metabolome_normalization.txt'))
HmtMetabolomicsDataClass <- R6::R6Class("HmtMetabolomicsDataClass",

  inherit = AbstractMetabolomicsDataClass,

  public = list(

    initialize = function(filename)
    {
      super$initialize(filename)
      print('HmtMetabolomicsDataClass')
      self$data <- read.csv(filename, header = TRUE, check.names = FALSE, row.names = 1)
      invisible(self)
    }
  )
)
