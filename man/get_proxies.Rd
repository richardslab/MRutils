% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proxy_search.R
\name{get_proxies}
\alias{get_proxies}
\title{Find proxy snps to the list provided.}
\usage{
get_proxies(
  rsids,
  token = Sys.getenv("LDLINK_TOKEN"),
  population,
  results_dir = NULL,
  skip_api = FALSE,
  r2_threshold = 0.9
)
}
\arguments{
\item{rsids}{a list of rsids (of snps) for which proxies are wanted.}

\item{token}{a token to the ldlink nih service. default will look for an environment variable
LDLINK_PROXY. If you need one go here: https://ldlink.nci.nih.gov/?tab=apiaccess}

\item{population}{specify within which population you want to find proxies. E.g. "CEU", "YRI", etc.}

\item{results_dir}{A subdirectory (a hash that depends on the population and r2_threshold values) of this directory will be used to
cache results}

\item{skip_api}{A boolean indicating whether to only use the cached results by skipping the API calls.}

\item{r2_threshold}{The R^2 threshold to use when returning results}
}
\description{
function will keep a cache of previously queried proxies in the results dir and only query
snps that are not cached. This includes rsids that errored. in that case, an empty file will be
stored in the cache directory and the snp will be skipped in a subsequent call. Note that network errors
may cause empty files to erroneously indicate that a snp has no rsID.
}
\details{
A token for the service is required, please get one at https://ldlink.nci.nih.gov/?tab=apiaccess
}
\examples{

\dontrun{ #because it needs a token.
# If one does `Sys.setenv(LDLINK_TOKEN=<Your token here>)` first, this will work

get_proxies("rs2495477", pop="CEU", results_dir="derived_data")  #returns one proxy
get_proxies("rs373780327", pop="CEU", results_dir="derived_data") #returns no proxies

get_proxies(c("rs2495477","rs373780327"), pop="CEU", results_dir="derived_data")  #returns one proxy
get_proxies(c(), pop="CEU", results_dir="derived_data")  #returns no proxy (empty list!)

# note that thanks to the cache, calling either of these a second time will result in an immediate
# return value

}

}
