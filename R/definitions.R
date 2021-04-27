

#' Validate a daaframe as a gwas
#' 
#' Validate that a dataframe contains values that are consistent with being a gwas.
#' This mean validating accoring to the validation rule in \code{\link{gwas_rules}}
#'
#' @param data input data, a dataframe
#'
#' @export
#'
assert_gwas <- function(data) {
  rsid <- CHR <- POS <- SE <- beta <- EAF <- NEA <- EA <- P <- NULL
  
  gwas_rules <- validate::validator(validate::field_format(rsid, "rs*"),
                                    validate::field_format(rsid, "rs[0-9]*", type = "regex"),
                                    CHR %in% valid_contigs,
                                    POS > 0,
                                    validate::field_format(EA, "[ACGT]", type = "regex"),
                                    validate::field_format(NEA, "[ACGT]", type = "regex"),
                                    validate::in_range(EAF, 0, 1),
                                    is.numeric(beta),
                                    is.numeric(SE),
                                    validate::in_range(SE, min = 0, Inf, strict = TRUE),
                                    validate::in_range(P, 0, 1, strict = TRUE)
  )
  
  validations <- summary(validate::confront(data, gwas_rules))
  assertthat::assert_that(!any(validations$error))
 
}


#' Title check that input is a vector of strings that look like rsids
#'
#' @param strings 
#'
#' @return whether the input is of type and contents consistent with being a co
#' @export
#'
#' @examples
#' assert_rsids(c("rs001101","rs00042")) # TRUE
#' assert_rsids(c("rs001101")) # TRUE
#' assert_rsids("rs001101") # TRUE
#' 
#' if(false) {
#'    assert_rsids(c("001101","rs00042")) ## error
#' }
#' 
assert_rsids <- function(strings) {
  assertthat::assert_that(typeof(data$rsid) == "character")
  assertthat::assert_that(all(grepl("^rs[0-9]*", strings)))
}


assert_probability <- function(p) {
  
  assertthat::assert_that(is.numeric(p))
  assertthat::assert_that(p >= 0 && p <= 1)
}

required_headers <- c("rsid", "CHR", "POS", "P", "beta", "EA", "NEA", "EAF", "SE")

valid_contigs <- c(1:22, "X", "Y")
valid_contigs <- c(valid_contigs, paste0("chr",valid_contigs))

valid_references <- c("hg18", "hg19", "hg38")

valid_alleles <- c("A", "C", "G", "T")

palindromic <- c("(A/T)", "(T/A)", "(C/G)", "(G/C)")

