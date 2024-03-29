% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsid_search.R
\name{merge_rsids_into_gwas}
\alias{merge_rsids_into_gwas}
\title{Method to merge back rsids from a subset of a gwas into a full gwas where they were missing}
\usage{
merge_rsids_into_gwas(gwas, rsids, validate = TRUE)
}
\arguments{
\item{gwas}{a dataframe containing the gwas with CHR POS NEA and EA}

\item{rsids}{another gwas which contained a subset of the rows in gwas, presumably with
some rsids updated.
#' @param validate a boolean indicating whether to validate the resulting gwas
(TRUE by default) use FALSE for debugging.}
}
\value{
original input gwas with available rsids replacing rsids where possible
}
\description{
Method to merge back rsids from a subset of a gwas into a full gwas where they were missing
}
\examples{

any(is.na(demo_data$rsid)) # TRUE
fixed_partial_gwas <- get_unknown_rsids_from_locus(demo_data)
any(is.na(fixed_partial_gwas$rsid)) # FALSE
nrow(fixed_partial_gwas) == nrow(demo_data) # FALSE
fixed_gwas <- merge_rsids_into_gwas(demo_data, fixed_partial_gwas)
nrow(fixed_gwas) == nrow(demo_data) # TRUE

still_demo_data <- merge_rsids_into_gwas(demo_data, subset(fixed_partial_gwas, FALSE)) # same as original
dplyr::all_equal(still_demo_data, demo_data) # TRUE

}
\keyword{internal}
