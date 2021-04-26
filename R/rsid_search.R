
#' Get rsid name from position
#' 
#' This method accesses the NCI's API and queries a position (CHR & POS) and variant (REF & ALT)
#' to obtain a rsid. If it fails, it will attempt to swap REF and ALT before giving up.
#'
#' @param chrom The name of the contig for the SNP
#' @param pos The (1-based) position of the SNP 
#' @param ref The reference allele
#' @param alt The alternate allele
#' @param assembly Which reference genome to use ("hg19" or "hg38")
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
get_rsid_from_position <- function(chrom, pos, ref, alt, assembly = "hg19") {
  retVal=tryCatch({
    baseURL1 <- "https://api.ncbi.nlm.nih.gov/variation/v0/vcf/{chrom}/{pos}/{ref}/{alt}/contextuals?assembly={assembly}"
    baseURL1_swapped <- "https://api.ncbi.nlm.nih.gov/variation/v0/vcf/{chrom}/{pos}/{alt}/{ref}/contextuals?assembly={assembly}"
    
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
    
    baseURL2 <- "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/{seq_id}:{pos}:{ref}:{alt}/rsids"
    baseURL2_swapped <- "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/{seq_id}:{pos}:{alt}:{ref}/rsids"
    
    id <- tryCatch({
      url <- glue::glue(baseURL2)
      Sys.sleep(1)
      paste0("rs", jsonlite::read_json(url)$data$rsids[[1]])
    },
    error = function(e){
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
    warning(paste("there was an error:",e))
    NULL
  }
  )
  as.character(retVal)
}


get_unknown_rsids_from_locus <- function(gwas, build = "hg19") {
  
  rsid <- CHR <- POS <- NEA <- EA <- NULL
  
  unknown_ids <- subset( x = gwas, 
                         subset = is.na(rsid), 
                         select = c(CHR, POS, NEA, EA)) %>% transform(build = build)
  
  ## this can take time and hits the API multiple times....
  withIds <- plyr::adply(unknown_ids, 1, function(x) {
    c(rsid = get_rsid_from_position(chrom = x$CHR, 
                                    pos = x$POS, 
                                    ref = x$NEA, 
                                    alt = x$EA, 
                                    assembly = x$build))
  }
  )
  withIds
}


merge_rsids_into_gwas <- function(gwas, rsids) {
  
  rsid <- CHR <- POS <- rsid.x <- rsid.y <- NULL
  
  gwas_with_ids <- plyr::mutate(gwas, rsid = as.character(rsid)) %>%
    merge(subset(rsids, select = c(CHR, POS, rsid)), all.x = T, by = c("CHR","POS")) %>%
    plyr::mutate(rsid = plyr::adply(do.call(cbind, list(rsid.y, rsid.x)), 1, function(x) dplyr::first(stats::na.omit(x)))$V1, 
           rsid.x = NULL, 
           rsid.y = NULL)
  gwas_with_ids
}


`%notin%` <- Negate(`%in%`)

