


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
#' to a data-frame that will be sanity-checked as being a gwas.
#' @param validate a boolean indicating whether to validate the resulting gwas 
#' (after applying the mapping function) (TRUE by default) use FALSE for debugging.
#' @param chr_action indicates whether to add a "chr" prefix (if not present) to the names of the contigs in the 
#' query SNPs, remove it (where present), leave as is, or try to do both (default). "both", "leave","remove","add".
#' @param pad the amound (in basepairs) to pad each provided locus prior to querying for the results.
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
#' system(glue::glue("tabix -s 1 -b 2 -e 2 -f '{partial_gwas}'"))
#' 
#' mapping_function <- function(x) {
#'    dplyr::mutate(x,
#'       CHR=as.character(X.CHR),
#'       beta = all_inv_var_meta_beta, 
#'       P = all_inv_var_meta_p, 
#'       EA = ALT, 
#'       NEA = REF, 
#'       EAF = all_meta_AF, 
#'       SE = all_inv_var_meta_sebeta) %>% subset(select = required_headers)
#' }
#' extract_snps_from_bgzip(partial_gwas, demo_data, mapping_function=mapping_function, chr_action="remove")
#'
extract_snps_from_bgzip <-
  function(
    outcome,
    snps,
    mapping_function = identity,
    comment_char = "",
    validate = TRUE,
    chr_action = "both",
    pad = 0
  ) {
    CHR <- POS <- . <- NULL 
    
    chr_action <-
      match.arg(chr_action, c("leave", "remove", "add", "both"))
    
    unique_proxies_raw <- unique(dplyr::select(snps, c("CHR", "POS")))
    if (chr_action == "add") {
      proxies_for_query <- unique_proxies_raw %>%
        dplyr::mutate(
          CHR = gsub(
            "^(?!chr)(?#non-present chr in start of string)",
            "chr",
            unique_proxies_raw$CHR,
            perl = TRUE
          )
        ) 
    } else if (chr_action == "remove") {
      proxies_for_query <- unique_proxies_raw %>%
        dplyr::mutate(CHR = gsub("^chr", "", CHR))
    } else if (chr_action == "both") {
      proxies_for_query <- rbind(
        unique_proxies_raw %>%
          dplyr::mutate(
            CHR = gsub("^(?!chr)", "chr",
                       unique_proxies_raw$CHR,
                       perl = TRUE)
          ),
        unique_proxies_raw %>%
          dplyr::mutate(CHR = gsub("^chr", "", CHR))
      )
      
    }
    proxies_for_query <- proxies_for_query %>% dplyr::select(c(CHR, POS))%>%dplyr::mutate(start=pmax(1,POS-pad), end=POS+pad)
    region <- with(proxies_for_query, glue::glue("{CHR}\t{start}\t{end}")) %>% paste(collapse ="\n") 
    # region <- with(proxies_for_query,n=1000, glue::glue("{CHR}\t{start}\t{end}")) %>% paste(collapse ="\n") 
    temp_region <- tempfile(pattern = "region__",
                            tmpdir = tempdir(),
                            fileext = ".txt")
    
    temp_region_file <- file(temp_region,"wt")
    cat(region, file=temp_region_file, fill = TRUE)
    close(temp_region_file)
    
    temp_output <- tempfile(pattern = "subsetted_gwas__",
                            tmpdir = tempdir(),
                            fileext = ".txt")

    # -f to force even if index is older
    # -T to use temp_region as "targets"
    cmd <- glue::glue("tabix -f -T '{temp_region}' '{outcome}' > '{temp_output}'")
    system(cmd)
    print(cmd)
    
    col_names <- utils::read.table(
      outcome,
      header = TRUE,
      nrows = 1,
      comment.char = comment_char,
    ) %>% names()
    
     outcome_data <-
      utils::read.table(temp_output, comment.char = "", col.names = col_names, as.is=TRUE) %>%
      mapping_function()
    if (validate) {
      assert_gwas(outcome_data)
    }
    outcome_data
  }



#' Function that is a mapping function for the example data in system.file("extdata", "COVID_partial_gwas_hg37.txt.gz", package = "MRutils")
#' Can be used when accessing that data for example with extract_snps_from_bgzip
#'
#' @name example_mapping_function
#'
#' @param x a dataframe which was read from a file
#' 
#' @return a dataframe that was mapped and should now be valid according to assert_gwas
#'
#' @export
#'
example_mapping_function <- function(x) {
  X.CHR<-all_inv_var_meta_beta<-all_inv_var_meta_p<-ALT<-REF<-all_meta_AF<-all_inv_var_meta_sebeta<-NULL

  dplyr::mutate(x,
               CHR = as.character(X.CHR),
               beta = all_inv_var_meta_beta,
               P = all_inv_var_meta_p,
               EA = ALT,
               NEA = REF,
               EAF = all_meta_AF,
               SE = all_inv_var_meta_sebeta) %>% 
    subset(select = required_headers)
}


