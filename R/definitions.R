

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
#' @name required_headers
#' @export
#' @keywords internal
valid_contigs <- c(1:22, "X", "Y")

#' List of valid names for a reference sequence/contig in this package
#'
#' This only allows for contigs with the "chr" prefix. For the list containing see
#' \code{\link{valid_contigs}}
#' @name required_headers
#' @export
#' @keywords internal
valid_contigs_with_chr <- paste0("chr", valid_contigs)

#' List of human reference builds that can be used to find rsids
#'
#' @name required_headers
#' @keywords internal
valid_references <- c("hg18", "hg19", "hg38")

#' Validate a dataframe as a gwas
#'
#' Validate that a dataframe contains values that are consistent with being a gwas.
#'
#' @param data input data, a dataframe
#' @export
#' @examples
#'
#'  assert_gwas(demo_data) # TRUE
#'
#'  \dontrun{
#'     broken_data <- demo_data # make copy
#'     broken_data$POS[1] <- 0 # Zero is not a valid value for POS
#'     assert_gwas(broken_data) # NOT OK. Will error.
#'  }
#'
#'
assert_gwas <- function(data, show_error=TRUE) {
  rsid <- CHR <- POS <- SE <- beta <- EAF <- NEA <- EA <- P <- NULL
  
  . <- field_format <- in_range <- NULL
  
  gwas_rules <-
    validate::validator(
      rsid_starts_rs = is.na(rsid) | field_format(rsid, "rs*"),
      rsid_has_numbers = is.na(rsid) |
        field_format(rsid, "^rs[0-9]*$", type = "regex"),
      chr_is_valid = CHR %in% valid_contigs |
        all(CHR %in% valid_contigs_with_chr),
      chr_is_valid_with_chr = CHR %in% valid_contigs_with_chr |
        all(CHR %in% valid_contigs),
      pos_is_positive = POS > 0,
      ea_is_dna = field_format(EA, "[ACGT]", type = "regex"),
      nea_is_dna = field_format(NEA, "[ACGT]", type = "regex"),
      eaf_is_prob = in_range(EAF, 0, 1),
      se_is_postive = in_range(SE, min = 0, Inf, strict = TRUE),
      p_is_prob = in_range(P, 0, 1, strict = TRUE),
      p_is_numeric = is.numeric(P),
      beta_is_numeric = is.numeric(beta),
      eaf_is_numeric = is.numeric(EAF),
      se_is_numeric = is.numeric(SE),
      pos_is_numeric = is.numeric(POS),
      required_headers_present = all(required_headers %in% names(.))
      
    )
  validate::voptions(gwas_rules, raise = 'all')
  
  
  val_sum <-
    validate::summary(validate::confront(as.data.frame(data), gwas_rules))
  if (any(val_sum$error) || any(val_sum$fails > 0)) {
    if(show_error) {
      methods::show(val_sum)
      methods::show(validate::violating(as.data.frame(data), gwas_rules[1:10]))
    }
    assertthat::assert_that(!any(val_sum$error))
    assertthat::assert_that(!any(val_sum$fails))
  }
  
  invisible(TRUE)
}


#' Check that input is a vector of strings that look like rsids
#'
#' @param strings input collection of strings
#'
#' @return whether the input is of type and contents consistent with being a
#' collection of rsids
#'
#' @export
#'
#' @examples
#' assert_rsids(c("rs001101","rs00042")) # TRUE
#' assert_rsids(c("rs001101")) # TRUE
#' assert_rsids("rs001101") # TRUE
#'
#' if (FALSE) {
#'    assert_rsids(c("001101","rs00042")) ## error
#' }
#'
assert_rsids <- function(strings) {
  assertthat::assert_that(typeof(strings) == "character")
  assertthat::assert_that(all(grepl("^rs[0-9]*", strings)))
}


#' Asserts that input is a numeric vector or scalar and that all the values are
#' between 0 and 1 inclusive
#'
#' @param p input scalar or vector to be verified
#'
#' @export
#' @examples
#'
#' assert_probability(0) # OK
#' assert_probability(1) # OK
#' assert_probability(c(0.0,.1,.04)) # OK
#'
#' if(FALSE) {
#'     assert_probability(2.4) # not OK, will error
#'     assert_probability(-0.1) # not OK, will error
#'     assert_probability(c(0.1,2.1,-0.2)) # not OK, will error
#' }
assert_probability <- function(p) {
  assertthat::assert_that(is.numeric(p))
  assertthat::assert_that(all(p >= 0 & p <= 1))
}
