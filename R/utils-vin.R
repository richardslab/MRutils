#' vin operator
#'
#' See \code{validate::\link[validate:vin]{\%vin\%}} for details.
#'
#' @name %vin%
#' @rdname vin
#' @keywords internal
#' @export
#' @importFrom validate %vin%
#' @usage x \%vin\% table
#' @param vector or NULL: the values to be matched
#' @param table vector or NULL: the values to be matched against.
#' @return A vector of the same length as x. TRUE if the element in x is 
#' definitely in table, FALSE if not and NA is unknown (e.g. if table contains NA)
NULL
