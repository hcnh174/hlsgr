AbstractBaseClass <- R6::R6Class("AbstractBaseClass",

  public = list(
    initialize = function(dds, outdir)
    {
      print('AbstractBaseClass')
      invisible(self)
    }
  )
)

AbstractDifferentialExpressionAnalysisClass <- R6::R6Class("AbstractDifferentialExpressionAnalysisClass",

  inherit = AbstractBaseClass,

  public = list(
    initialize = function()
    {
      super$initialize()
      print('AbstractDifferentialExpressionAnalysisClass')
      invisible(self)
    }
  )
)




