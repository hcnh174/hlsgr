#' Reassure A User That They Are Doing Well
#'
#' Look at that, they're building a documented function already!
#'
#' @param me the user to reassure
#'
#' @return NULL (invisibly). Used for the side-effect of generating a
#' \code{message}.
#' @export
#'
#' @examples
#' makeMeSmile("You")
makeMeSmile <- function(me = "I") {
  message(me, " can develop packages now!")
}
