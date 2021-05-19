

#' Get rsid name from position
#'
#' This method combines the cached and API version of get_rsid 
#' 
#' @param chrom The name of the contig for the SNP
#' @param pos The (1-based) position of the SNP
#' @param ref The reference allele (only SNPs are supported)
#' @param alt The alternate allele (only SNPs are supported)
#' @param assembly Which reference genome to use ("hg18", "hg19", or "hg38")
#' @param cache_file Which file to use for caching results. if NULL will not use cache
#' @param update_cache A boolean indicating whether to update the cache results.
#'
#' @return an rsid identifier of the position provided
#' 
#' @keywords internal
#' @export
#'
#' @examples
#'
#'get_rsid_from_position("9", 125711603,	"C",	"A") # "rs10760259"
#'
#'
#'
get_rsid <-
  function(chrom,
           pos,
           ref,
           alt,
           assembly = valid_references,
           cache_file = NULL,
           update_cache = TRUE) {
    #try cache?
    assembly <- match.arg(assembly)
    try_cache <- is.null(cache_file) || file.exists(cache_file)
    update_cache <- try_cache && update_cache
    
    rsid = NULL
    if (try_cache) {
      print(glue::glue("getting cached rsid for {chrom}:{pos}"))
      rsid = get_cached_rsid_from_position(chrom, pos, ref, alt, assembly = "hg19")
    }
    
    if (is.null(rsid)) {
      print(glue::glue("getting fresh rsid for {chrom}:{pos}"))
      rsid = get_rsid_from_position(chrom, pos, ref, alt, assembly = "hg19")
    } else {
      update_cache = FALSE
    }
    if (!is.null(rsid) && update_cache) {
      print(glue::glue("writing rsid to cache {rsid}"))
      Sys.sleep(1)
      put_rsid_into_cache(rsid, chrom, pos, ref, alt, assembly, cache_file)
    }
    
    rsid
  }


#' Get rsid name from position
#'
#' This method accesses the NCI's API and queries a position (CHR & POS) and variant (REF & ALT)
#' to obtain a rsid. If it fails, it will attempt to swap REF and ALT before giving up.
#'
#' @param chrom The name of the contig for the SNP
#' @param pos The (1-based) position of the SNP
#' @param ref The reference allele (only SNPs are supported)
#' @param alt The alternate allele (only SNPs are supported)
#' @param assembly Which reference genome to use ("hg18", "hg19", or "hg38")
#'
#' @return an rsid identifier of the position provided
#' @export
#'
#' @examples
#'
#'get_rsid_from_position("9", 125711603,	"C",	"A") # "rs10760259"
#'
#'
#'
get_rsid_from_position <-
  function(chrom, pos, ref, alt, assembly = valid_references) {
    assembly <- match.arg(assembly)
    assertthat::assert_that(assembly %in% valid_references)
    assertthat::assert_that(chrom %in% valid_contigs |
                              chrom %in% valid_contigs_with_chr)
    assertthat::assert_that(ref %in% valid_alleles)
    assertthat::assert_that(alt %in% valid_alleles)
    assertthat::assert_that(pos > 0)
    
    retVal <- tryCatch({
      baseURL1 <-
        "https://api.ncbi.nlm.nih.gov/variation/v0/vcf/{chrom}/{pos}/{ref}/{alt}/contextuals?assembly={assembly}"
      baseURL1_swapped <-
        "https://api.ncbi.nlm.nih.gov/variation/v0/vcf/{chrom}/{pos}/{alt}/{ref}/contextuals?assembly={assembly}"
      
      f <- tryCatch({
        url <- glue::glue(baseURL1)
        Sys.sleep(1)
        jsonlite::read_json(url)$data$spdis[[1]]
      },
      error = function(e) {
        warning("There was an error (1):")
        warning(e)
        warning("Trying to swap ref and alt")
        Sys.sleep(1)
        jsonlite::read_json(glue::glue(baseURL1_swapped))$data$spdis[[1]]
      })
      
      #nolint (unused variable)
      pos <- f$position
      #nolint (unused variable)
      seq_id <- f$seq_id
      
      baseURL2 <-
        "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/{seq_id}:{pos}:{ref}:{alt}/rsids"
      baseURL2_swapped <-
        "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/{seq_id}:{pos}:{alt}:{ref}/rsids"
      
      id <- tryCatch({
        url <- glue::glue(baseURL2)
        Sys.sleep(1)
        paste0("rs", jsonlite::read_json(url)$data$rsids[[1]])
      },
      error = function(e) {
        warning("There was an error (2):")
        warning(e)
        warning("Trying to swap ref and alt")
        url <- glue::glue(baseURL2_swapped)
        Sys.sleep(1)
        id <- jsonlite::read_json(url)$data$rsids[[1]]
        glue::glue("rs{id}")
      })
    },
    error = function(e) {
      warning(paste("there was an error:", e))
      NULL
    })
    retVal
  }

col_types <- readr::cols(
  RSID = readr::col_character(),
  CHR = readr::col_double(),
  POS = readr::col_double(),
  REF = readr::col_character(),
  ALT = readr::col_character()
)


#' Get rsid name from position
#'
#' This method accesses the NCI's API and queries a position (CHR & POS) and variant (REF & ALT)
#' to obtain a rsid. If it fails, it will attempt to swap REF and ALT before giving up.
#'
#' @param chrom The name of the contig for the SNP
#' @param pos The (1-based) position of the SNP
#' @param ref The reference allele (only SNPs are supported)
#' @param alt The alternate allele (only SNPs are supported)
#' @param assembly Which reference genome to use ("hg18", "hg19", or "hg38")
#'
#' @return an rsid identifier of the position provided
#' 
#' @keywords internal
#' @export
#'
get_cached_rsid_from_position <-
  function(chrom,
           pos,
           ref,
           alt,
           assembly = valid_references,
           cache = NULL) {
    assembly <- match.arg(assembly)
    if (is.null(cache)) {
      cache_file = here::here(glue::glue("cache/{assembly}.dbSNP"))
    }
    
    CHR <- POS <- REF <- ALT <- NULL
    
    if (file.exists(cache_file)) {
      cached = readr::read_tsv(cache_file, col_types = col_types, progress = FALSE) %>% 
        subset(CHR == chrom &
                 POS == pos &
                 (REF == ref &
                    ALT == alt |
                    REF == alt &
                    ALT == ref))
      if (nrow(cached) == 1)
        return(cached$RSID[1])
      if (nrow(cached) > 1)
        stop(glue::glue("found more than one row corresponding to query: {cached}"))
    } else {
      cat(glue::glue("cache file {cache_file} doesn't exist."))
    }
    return(NULL)
  }



put_rsid_into_cache <-
  function(rsid,
           chrom,
           pos,
           ref,
           alt,
           assembly = valid_references,
           cache = NULL) {
    assembly <- match.arg(assembly)
    if (is.null(cache)) {
      cache_file = here::here(glue::glue("cache/{assembly}.dbSNP"))
    }
    
    RSID <- CHR <- POS <- REF <- ALT <- NULL
    
     
    new_row = data.frame(
      RSID = rsid,
      CHR = chrom,
      POS = pos,
      REF = ref,
      ALT = alt
    )
    
    if (!file.exists(dirname(cache_file))) {
      dir.create(dirname(cache_file), recursive = TRUE)
    }
    
    if (!file.exists(cache_file)) {
      
      readr::write_tsv(new_row, cache_file)
      return(FALSE)
    } else {
      cached = get_cached_rsid_from_position(chrom, pos, ref, alt)
      
      if (is.null(cached)) {
        readr::write_tsv(new_row, cache_file, append = TRUE)
        return(FALSE)
      }
      if (cached != rsid)
        stop(glue::glue("found row with same coordinates but different id: {cached}"))
    }
    readr::write_tsv(new_row, cache_file, append = TRUE)
    return(TRUE)
  }

#' Method to obtain rsids for records in a gwas where they are missing
#'
#' The method looks at all the rows for which rsid is NA and gets their rsid by
#' calling get_unknown_rsids_from_locus and merging the resulting rsids into the gwas
#' using merge_rsids_into_gwas
#'
#' @name fill_gwas_unknown_rsids
#' @param gwas a dataframe containing the gwas with CHR POS NEA and EA
#' @param assembly one of "hg18", "hg19" (default), "hg38".
#'
#' @return original input gwas with available rsids replacing rsids where possible
#' @export
#'
#' @examples
#'
#' any(is.na(demo_data$rsid)) # TRUE
#' fixed_gwas = fill_gwas_unknown_rsids(demo_data)
#' any(is.na(fixed_gwas$rsid)) # FALSE
#' nrow(fixed_gwas) == nrow(demo_data) # TRUE
#' 
#'
fill_gwas_unknown_rsids <- function(gwas, assembly = valid_references) {
  assembly <- match.arg(assembly)
  assert_gwas(gwas)
  withIds <- get_unknown_rsids_from_locus(gwas, assembly)
  merge_rsids_into_gwas(gwas, withIds)
}


#' Method to obtain rsids for records in a gwas where they are missing
#'
#' The method looks at all the rows for which rsid is NA and gets their rsid by
#' calling merge_rsids_into_gwas and merging the resulting rsids into the gwas
#'
#'
#' @param gwas a dataframe containing the gwas with CHR POS NEA and EA
#' @param assembly one of "hg18", "hg19" (default), "hg38".
#'
#' @return subset of input gwas with available rsids replacing rsids where they were originally missing. Does
#' not include all the columns in the original gwas, only CHR, POS, NEA, EA, rsid
#' @export
#' @keywords internal
#'
#'
#' @examples
#'
#' any(is.na(demo_data$rsid)) # TRUE
#' partial_gwas = get_unknown_rsids_from_locus(demo_data)
#' any(is.na(partial_gwas$rsid)) # FALSE
#' nrow(partial_gwas) == nrow(demo_data) # false
#' 
#'
get_unknown_rsids_from_locus <- function(gwas, assembly = valid_references) {
  assembly = match.arg(assembly)
  assert_gwas(gwas)
  assertthat::assert_that(assembly %in% valid_references)
  
  rsid <- CHR <- POS <- NEA <- EA <- NULL
  
  unknown_ids <- subset(
    x = gwas,
    subset = is.na(rsid),
    select = c(CHR, POS, NEA, EA)
  ) %>%plyr::mutate(assembly=assembly)
  
  ## this can take time and hits the API multiple times....
  withIds <- plyr::adply(unknown_ids, 1, function(x) {
    c(rsid = get_rsid(
      chrom = x$CHR,
      pos = x$POS,
      ref = x$NEA,
      alt = x$EA,
      assembly = x$assembly
    ))
  }) %>% plyr::mutate(assembly=NULL)
  
  withIds
}



#' Method to merge back rsids from a subset of a gwas into a full gwas where they were missing
#'
#'
#' @param gwas a dataframe containing the gwas with CHR POS NEA and EA
#' @param rsids another gwas which contained a subset of the rows in gwas, presumably with 
#' some rsids updated. 
#'
#' @return original input gwas with available rsids replacing rsids where possible
#' @export
#'
#' @keyword internal
#' @examples
#'
#' any(is.na(demo_data$rsid)) # TRUE
#' fixed_partial_gwas <- get_unknown_rsids_from_locus(demo_data)
#' any(is.na(fixed_partial_gwas$rsid)) # FALSE
#' nrow(fixed_partial_gwas) == nrow(demo_data) # FALSE
#' fixed_gwas <- merge_rsids_into_gwas(demo_data, fixed_partial_gwas)
#' nrow(fixed_gwas) == nrow(demo_data) # TRUE
#' 
merge_rsids_into_gwas <- function(gwas, rsids) {
  assert_gwas(gwas)
  assert_rsids(rsids)
  
  rsid <- CHR <- POS <- rsid.x <- rsid.y <- NULL
  
  gwas_with_ids <- plyr::mutate(gwas, rsid = as.character(rsid)) %>%
    merge(subset(rsids, select = c(CHR, POS, rsid)),
          all.x = T,
          by = c("CHR", "POS")) %>%
    plyr::mutate(
      rsid = plyr::adply(do.call(cbind, list(rsid.y, rsid.x)), 1, function(x)
        dplyr::first(stats::na.omit(x)))$V1,
      rsid.x = NULL,
      rsid.y = NULL
    )
  gwas_with_ids
}

#' Not in operator
#'
#' Negates the result of %in%
#'
#' @name %notin%
#' @rdname not_in
#' @keywords internal
#' @usage lhs \%notin\% rhs
#' @param lhs vector or NULL: the values to be matched. Long vectors are supported.
#' @param rhs vector or NULL: the values to be matched. Long vectors are not supported.
#' @return A logical vector, indicating if a match was NOT located for each element of x: thus the values are TRUE or FALSE and never NA.
`%notin%` <- function(lhs, rhs) {
  match(lhs, rhs, nomatch = 0) == 0
}
