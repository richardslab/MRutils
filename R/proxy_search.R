


#' Find proxy snps to the list provided. 
#' 
#' function will keep a cache of previously obtained proxies in the results dir and only query 
#' snps that are not cached. A token for the service is required, please get one at https://ldlink.nci.nih.gov/?tab=apiaccess
#'
#'
#' @name get_proxies
#'
#' @param rsid a list of rsids (of snps) for which proxies are wanted.  
#' @param token a token to the ldlink nih service. default will look for an environment variable 
#' LDLINK_PROXY. If you need one go here: https://ldlink.nci.nih.gov/?tab=apiaccess
#' @param population specify within which population you want to find proxies. E.g. "CEU", "YRI", etc. go 
#' @param results_dir A subdirectory (a hash that depends on the population and r2_threshold values) of this directory will be used to 
#' cache results
#' @param skip_api A boolean indicating whether to only use the cached results by skipping the API calls.
#' @param r2_threshold The R^2 threshold to use when returning results 
#'
#' @export
#'
#'
get_proxies <- function(rsids, token = Sys.getenv("LDLINK_TOKEN"), population, results_dir, skip_api = FALSE, r2_threshold = 0.9) {
  config=list(pop=population, r2_thresh=r2_threshold)
  hashed_subdir=digest::digest(config)
  hashed_subdir_full=glue::glue("{results_dir}/{hashed_subdir}")
  
  if(!dir.exists(hashed_subdir_full)) {
    dir.create(hashed_subdir_full)
    xfun::in_dir(hashed_subdir_full, {yaml::write_yaml(config,"config.yaml")} )
  }
  
  xfun::in_dir(hashed_subdir_full, {
    
    raw_rs_files=dir(".", "^rs[0-9]*\\.txt")
    
    rsFiles <- Filter(function(x) file.info(x)$size > 1000, raw_rs_files)

    missing_rsids <- setdiff(rsids, gsub("\\.txt$","", raw_rs_files)) 
    
    # this takes time and hits the LDlink API.
    if (length(missing_rsids)>0) {
      if( !skip_api) {
        LDlinkR::LDproxy_batch(missing_rsids, pop = population, r2d = "r2", append = F, token = token)
        rsFiles <- Filter(function(x) file.info(x)$size > 1000, dir(".", "^rs[0-9]*\\.txt"))
      } else {
        warning(glue::glue("Not calling the LDproxy api, but there are some missing rsids in the results dir:{missing_rsids}"))
      }
    }
    
    
    R2 <- Locus <- Coord <- Alleles <- RS_Number <- rsid <- query_rsid <- NULL
    # only read snps  
    proxies <- plyr::ldply(rsFiles, function(x) dplyr::mutate(query_rsid = gsub("\\.txt$", "", x), utils::read.table(x, sep = "\t")) %>% 
                       subset(R2 >= 0 & grepl("([ACGT]/[ACGT])", x = Alleles))) %>% 
      dplyr::mutate(rsid = RS_Number, Locus = Coord) %>% 
      subset(query_rsid %in% rsids) %>%
      dplyr::select(c(-"Distance", -"Dprime", -"RegulomeDB", -"Function", -"RS_Number", -"Coord")) %>% 
      dplyr::mutate(CHR = gsub(pattern = "^chr", replacement = "", do.call(rbind, strsplit(as.character(Locus), split = ':', fixed = TRUE))[,1]),
             POS = as.integer(do.call(rbind, strsplit(as.character(Locus), split = ':', fixed = TRUE))[,2])) %>% 
      subset(R2 >= r2_threshold) %>% mutate_cond(condition = rsid==".", rsid=NA)
  })
  assert_proxies(proxies)
  proxies
}



#' @param rsids_and_chr is a dataframe that contains (at least) the following two columns: CHR, rsid
#' 

#' @param token is an access token to the nci API. 
#' if you don't have a token go here: https://ldlink.nci.nih.gov/?tab=apiaccess
#' @return
#' @export
#'
#' @examples
#' 
#' 
#' 
get_LD_pairs <- function(rsids_and_chr, population, token) { 
  assert_rsids()
  chr_validator
  CHR <- NULL
  # call LD matrix per each chromosome
  plyr::ddply(rsids_and_chr, .progress = "text", .variables = CHR, function(x) {
    if (nrow(x) > 1) {
      LDlinkR::LDmatrix(x$rsid, pop = population, token, r2d = "r2") %>%  reshape2::melt(id.vars = "RS_number")
    } else {
      temp <- data.frame(RS_number = as.character(x$rsid[1]), X = 1)
      names(temp)[2] <- as.character(x$rsid[1])
      temp %>% reshape2::melt(id.vars = "RS_number")
    }
  })
}



#' method to remove palindromicly ambiguous snps and select for the smallest pvalue within each "query_rsid":
#'
#' @param choose_best_proxies the dataframe containing results from LDprxy merged with an outcome gwas.
#' The Alleles column is assumed to have been already "fixed" according to the "Correlated_Alleles" column.
#' @param near_half_threshold The distance from AEF=0.5 that is to be considered "uimbiguous" when a palindromic SNP
#' has it.
#'
#' @return The original dataframe, with palindromic SNPs removed and within each query_rsid only the one with the smallest 
#' P value retained. In the case of a tie, the one with the smallest distance to the query snp (located in POS.exp) will be 
#' chosen. In the extremely rare case of a further tie, it will be broken randomly. 
#' @export
#'
#' @examples
choose_best_proxies <- function(proxies_and_more, near_half_threshold = 0.08) {
  proxies_and_more$POS=proxies_and_more$POS.exp
  assert_probabilities(proxies_and_more)
  assert_gwas(proxies_and_more)
  
  EAF <- POS.exp <- POS.proxy <- Alleles <- query_rsid <- maf_near_half <- R2 <- distance_to_query <- NULL
  
  proxies_and_more %>%
    dplyr::mutate(
      maf_near_half = abs(EAF - 0.5) <= near_half_threshold,
      distance_to_query = abs(POS.exp - POS.proxy)
    ) %>% 
    subset(!(Alleles %in% palindromic & maf_near_half)) %>%
    dplyr::group_by(query_rsid) %>%
    dplyr::slice_max(order_by = R2,
                     n = 1,
                     with_ties = TRUE) %>%
    dplyr::slice_min(order_by = distance_to_query,
                     n = 1,
                     with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(maf_near_half = NULL, distance_to_query = NULL, POS = NULL)
}


#' Prunes a collection of SNPs provided in a dataframe. 
#' The method is as follows: 
#'
#' @name prine_snps
#'
#' @param rsids_and_chr 
#' @param ld_pairs 
#' @param r2_threshold 
#'
#' @return a subset of the original dataframe, pruned according to some method....
#' @export
#'
#' @examples
prune_snps <- function(rsids_and_chr, ld_pairs, r2_threshold = 0.05) {
  
  RS_number <- rsid <- value <- variable <- NULL
  
  LDpairs_culled <- ld_pairs
  pairs <- subset(LDpairs_culled, value >= r2_threshold & as.character(RS_number) != as.character(variable))
  removed <- c()
  while (nrow(pairs) > 0) {
    ##TODO: make sure to keep snps with smaller P value.....
    to_remove <- names(which.max(table(c(as.character(pairs$RS_number),as.character(pairs$variable)))))
    # to_remove <- pairs[which.max(pairs$value), "variable"] %>% as.character()
    removed <- c(removed, to_remove)
    
    LDpairs_culled <- subset(LDpairs_culled, variable != to_remove & RS_number != to_remove)
    pairs <- subset(LDpairs_culled, value >= r2_threshold & as.character(RS_number) != as.character(variable))
  }
  
  list(rsids = subset(rsids_and_chr, rsid %notin% removed), 
       removed_rsid = removed)
}


#' Find proxy snps to the list provided. 
#' 
#' function will keep a cache of previously obtained proxies in the results dir and only query 
#' snps that are not cached. A token for the service is required, please get one at https://ldlink.nci.nih.gov/?tab=apiaccess
#'
#'
#' @name get_proxies
#'
#' @param rsid a list of rsids (of snps) for which proxies are wanted.  
#' @param token a token to the ldlink nih service. default will look for an environment variable 
#' LDLINK_PROXY. If you need one go here: https://ldlink.nci.nih.gov/?tab=apiaccess
#' @param population specify within which population you want to find proxies. E.g. "CEU", "YRI", etc. go 
#' @param results_dir A subdirectory (a hash that depends on the population and r2_threshold values) of this directory will be used to 
#' cache results
#' @param skip_api A boolean indicating whether to only use the cached results by skipping the API calls.
#' @param r2_threshold The R^2 threshold to use when returning results 
#'
#' @export
#'
#'
replace_alleles <- function(alleles, replacement_string) {
  . <- NULL
  allele_matrix <- strsplit(as.character(replacement_string), ",", fixed = T)[[1]] %>% 
  {plyr::llply(strsplit(x = ., split = "=", fixed = T))} %>% 
  {purrr::transpose(.l = .)} %>% 
    do.call(what = rbind)
  
  replace <- unlist(allele_matrix[2,])
  pattern <- unlist(allele_matrix[1,])
  
  a <- alleles %>% sapply(function(x) mgsub::mgsub(as.character(x), pattern = pattern, replacement = replace))
  a
}
