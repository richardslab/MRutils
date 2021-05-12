

#' List of genotypes that are palindromic
#'
#' @name palindromic
#' @keywords internal
palindromic <- c("(A/T)", "(T/A)", "(C/G)", "(G/C)")

#' List of strings that are valid alleles
#'
#' @name valid_alleles
#' @keywords internal
valid_alleles <- c("A", "C", "G", "T")

#' List of column headers that are required to be present in a GWAS (for this package)
#'
#' @name required_headers
#' @export
#' @keywords internal

required_headers <-
  c("rsid", "CHR", "POS", "P", "beta", "EA", "NEA", "EAF", "SE")

#' List of valid names for a reference sequence/contig in this package
#'
#' This only allows for contigs without the "chr" prefix. for the list containing chr see
#' \code{\link{valid_contigs_with_chr}}
#' @name valid_contigs
#' @export
#' @keywords internal
valid_contigs <- c(1:22, "X", "Y")

#' List of valid names for a reference sequence/contig in this package
#'
#' This only allows for contigs with the "chr" prefix. For the list containing see
#' \code{\link{valid_contigs}}
#' @name valid_contigs_with_chr
#' @export
#' @keywords internal
valid_contigs_with_chr <- paste0("chr", valid_contigs)

#' locus information validator
#' @name locus_validator
#' @keywords internal
#' 
locus_validator <-
  validate::validator(
    chr_is_valid = CHR %in% valid_contigs | all(CHR %in% valid_contigs_with_chr),
    chr_is_valid_with_chr = CHR %in% valid_contigs_with_chr |
    all(CHR %in% valid_contigs),
    pos_is_positive = POS > 0
    )

#' rsid information validator
#' @name rsid_validator
#' @keywords internal
rsid_validator <- validate::validator(
  rsid_starts_rs = is.na(rsid) | field_format(rsid, "rs*"),
  rsid_has_numbers = is.na(rsid) |field_format(rsid, "^rs[0-9]*$", type = "regex")
)


#' gwas column type validator
#' @name gwas_types_validator
#' @keywords internal
gwas_types_validator <- validate::validator(
  p_is_numeric = is.numeric(P),
  beta_is_numeric = is.numeric(beta),
  eaf_is_numeric = is.numeric(EAF),
  se_is_numeric = is.numeric(SE),
  pos_is_numeric = is.numeric(POS),
  chr_is_chr = is.character(CHR),
  required_headers_present = all(required_headers %in% names(.))
)

#' allele information validator
#' @name allele_validator
#' @keywords internal
allele_validator <- validate::validator(
  ea_is_dna = field_format(EA, "[ACGT]", type = "regex"),
  nea_is_dna = field_format(NEA, "[ACGT]", type = "regex"),
  eaf_is_prob = in_range(EAF, 0, 1)
)

proxy_allele_validator <- validate::validator(
  Alleles = field_format(Alleles, "([ACGT]/[ACGT])", type = "regex"),
  Correlated_Alleles = field_format(Correlated_Alleles, "[ACGT]=[ACGT],[ACGT]=[ACGT]", type = "regex")
)



#' stats information validator
#' @name stats_validator
#' @keywords internal
stats_validator <- validate::validator(
  se_is_postive = in_range(SE, min = 0, Inf, strict = TRUE),
  p_is_prob = in_range(P, 0, 1, strict = TRUE)
)

#' stats information validator
#' @name proxy_stats_validator
#' @keywords internal
proxy_stats_validator <- validate::validator(
  MAF_is_prob = in_range(MAF, min = 0, 1, strict = FALSE),
  R2_is_prob = in_range(Rs, 0, 1, strict = FALSE)
)


proxy_validator <- rsid_validator + locus_validator + proxy_stats_validator

#' locus information validator
#' @name gwas_validator
#' @keywords internal
gwas_validator <-
  rsid_validator +
  locus_validator +
  allele_validator +
  stats_validator

#' Check that input is a dataframe that validates accordig to the validator and optionally
#' show which data doesn't validate
#'
#' @param data a dataframe which will be validated against the validator
#' @param validator a validator against which to validate the data
#' @param show_error if data does _not_ validate, whether to show the reasons
#' 

assert_valid_data <- function(data, validator, show_error = c("all","none","summ") ) {
  show_error <- match.arg(show_error)
  val_sum <-
    validate::summary(validate::confront(data, validator))
  
  if (any(val_sum$error) || any(val_sum$fails > 0)) {
    if (show_error == "all" || show_error == "summ") {
      methods::show(subset(val_sum, fails!=0 | error))
    }
    if (show_error == "all"){
      methods::show(validate::violating(as.data.frame(data), validator))
      methods::show(validate::errors(as.data.frame(data), validator))
    }
    assertthat::assert_that(FALSE, "There's a problem with the data")
  }
  
  invisible(TRUE)
}



#' List of human reference builds that can be used to find rsids
#'
#' @name valid_references
#' @keywords internal
valid_references <- c("hg18", "hg19", "hg38")

#' Validate a dataframe as a gwas
#'
#' Validate that a dataframe contains values that are consistent with being a gwas.
#'
#' @param data input data, a dataframe
#' @param show_error if data does _not_ validate, whether to show the reasons
#' @export
#' @examples
#'
#'  assert_gwas(demo_data) # TRUE
#'  
#'  #TODO: add demo data with Rmpfr
#'
#'  \dontrun{
#'     broken_data <- demo_data # make copy
#'     broken_data$POS[1] <- 0 # Zero is not a valid value for POS
#'     assert_gwas(broken_data, show_error=TRUE) # NOT OK. Will error.
#'  }
#'
#'
assert_gwas <- function(data, show_error = c("all", "none", "summ")) {
  show_error <- match.arg(show_error)
  assert_valid_data(data, gwas_validator, show_error)
  assert_valid_data(data, gwas_types_validator, if(show_error=="all") "summ" else show_error)
}

#' Check that input is a vector of strings that look like rsids
#'
#' @param data a dataframe that has a column rsid which will be validated
#' @param show_error if data does _not_ validate, whether to show the reasons
#' 
#'
#' @return TRUE if valid, will throw a validation error is invalid
#' 
#'
#' @export
#'
#' @examples
#' assert_rsids(data.frame(rsid=c("rs001101","rs00042"))) # TRUE
#' assert_rsids(data.frame(rsid=c("rs001101"))) # TRUE
#' assert_rsids(data.frame(rsid="rs001101")) # TRUE
#'
#' if (FALSE) {
#'    assert_rsids(data.frame(rsid=c("001101","rs00042"))) ## error
#' }
#'
#'
assert_rsids <- curry::partial(assert_valid_data, list(validator=rsid_validator))


#' Check that input is a vector of strings that look like rsids
#'
#' @param data a dataframe that has a column rsid which will be validated
#' @param show_error if data does _not_ validate, whether to show the reasons
#' 
#'
#' @return TRUE if valid, will throw a validation error is invalid
#' 
#'
#' @export
#'
#' @examples
#' assert_rsids(data.frame(rsid=c("rs001101","rs00042"))) # TRUE
#' assert_rsids(data.frame(rsid=c("rs001101"))) # TRUE
#' assert_rsids(data.frame(rsid="rs001101")) # TRUE
#'
#' if (FALSE) {
#'    assert_rsids(data.frame(rsid=c("001101","rs00042"))) ## error
#' }
#'
#'
assert_probabilities <- curry::partial(assert_valid_data, list(validator=rsid_validator))

assert_proxies <- curry::partial(assert_valid_data, list(validator=proxy_validator))

assert_probability <- function(p) {
  assertthat::assert_that(p <= 1)
  assertthat::assert_that(0 <= p)
}


