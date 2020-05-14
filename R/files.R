#' Title
#'
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
checkFileExists <- function(filename)
{
  if (!file.exists(filename))
    R.oo::throw(concat('file does not exist: ',filename))
}
#checkFileExists('test.txt')
