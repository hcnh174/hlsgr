#' Concatenate strings similar to paste0
#'
#' @param ... strings to concatenate
#' @param sep separator, defaults to empty string
#' @param collapse defaults to NULL
#'
#' @return concatenated string
#' @export
#'
#' @examples
#' concat('abc','def','ghi')
concat <- function(..., sep='', collapse=NULL)
{
  strings<-list(...)
  #catch NULLs, NAs
  if (all(unlist(plyr::llply(strings,length))>0) &&	all(!is.na(unlist(strings))))
  {
    do.call("paste", c(strings, list(sep = sep, collapse = collapse)))
  }
  else
  {
    NULL
  }
}

#' Trims whitespace from string
#'
#' @param str input string to trim
#'
#' @return trimmed string
#' @export
#'
#' @examples
#' trim(' 	 abc def 	 ')
trim <- function(str)
{
  return(gsub("^\\s+|\\s+$", "", str))
}
