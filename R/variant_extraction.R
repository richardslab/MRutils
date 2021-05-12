


#' Extract a subset of a block-gzipped (and indexed) file by CHR and POS
#'
#' it is expected that snps has two columns 'CHR' and 'POS'
#' indicating the chromosome and the  (1-based) position of the snp in question.
#' It is assumed that `tabix` is installed on the system and directly available
#' for use by a system call.
#'
#' @param outcome gzipped file containing genomic positions and data. should be
#' already indexed by tabix
#' @param snps a dataframe containing at least the columns  'CHR' and 'POS'.
#' These position will be queried in the outcome
#' @param comment_char The character that is the comment character in the file.
#' @param mapping_function A function that will map the gwas as it is read from disk
#' to a dataframe that will be sanity-checked as being a gwas.
#'
#' @return a gwas containing the subset of the outcome file that was requested.
#'  The first non-commented line will be assumed to be the header line. It will be
#'  retained and used as the column names for the resulting sub-gwas.
#'
#' @export
#'
#' @examples
#'
#' partial_gwas <- system.file("extdata", "COVID_partial_gwas_hg37.txt.gz", package = "MRutils")
#' system(glue::glue("tabix -s 1 -b 2 -e 2 -f {partial_gwas}"))
#' 
#' mapping_function <- function(x) {
#'    dplyr::mutate(x,
#'       CHR=X.CHR,
#'       beta = all_inv_var_meta_beta, 
#'       P = all_inv_var_meta_p, 
#'       EA = ALT, 
#'       NEA = REF, 
#'       EAF = all_meta_AF, 
#'       SE = all_inv_var_meta_sebeta) %>% subset(select = required_headers)
#' }
#' extract_snps_from_bgzip(partial_gwas, demo_data, mapping_function=mapping_function)
#'
extract_snps_from_bgzip <-
  function(outcome,
           snps,
           mapping_function = identity,
           comment_char = "",
           validate = TRUE) {
    . <- NULL
    
    region <- snps %>%
      plyr::alply(1, function(x)
        with(x, glue::glue("{CHR}:{POS}-{POS}"))) %>% {
          do.call(paste, .)
        }
    
    col_names <- utils::read.table(
      outcome,
      header = TRUE,
      nrows = 1,
      comment.char = comment_char,
    ) %>% names()
    
    cat("column names", col_names)
    
    temp_output <- tempfile(pattern = "subsetted_exposure__",
                            tmpdir = tempdir(),
                            fileext = ".txt")
    
    cmd <- glue::glue("tabix -f '{outcome}' {region} > {temp_output}")
    system(cmd)
    outcome_data <-
      utils::read.table(temp_output, col.names = col_names) %>%
      mapping_function
    if (validate) {
      assert_gwas(outcome_data)
    }
    outcome_data
  }

