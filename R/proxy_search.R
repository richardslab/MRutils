

get_proxies <- function(rsids, token, population, results_dir, skip_api = FALSE, r2_threshold = 0.9) {
  
  xfun::in_dir(results_dir, {
    # this takes time and hits the LDlink API.
    if (!skip_api) {
      LDlinkR::LDproxy_batch(rsids, pop = population, r2d = "r2", append = F, token = token)
    }
    
    ## the filter is here because when LDproxy_batch failes, it creates a small file with the error message...
    rsFiles <- Filter(function(x) file.info(x)$size > 1000, 
                      dir(".", "^rs[0-9]*\\.txt"))
    
    
    R2 <- Locus <- Coord <- Alleles <- RS_Number <- NULL
    # only read snps  
    proxies <- plyr::ldply(rsFiles, function(x) dplyr::mutate(query_rsid = gsub("\\.txt$", "", x), utils::read.table(x, sep = "\t")) %>% 
                       subset(R2 >= 0 & grepl("([ACGT]/[ACGT])", x = Alleles))) %>% 
      dplyr::mutate(rsid = RS_Number, Locus = Coord) %>% 
      dplyr::select(c(-"Distance", -"Dprime", -"RegulomeDB", -"Function", -"RS_Number", -"Coord")) %>% 
      dplyr::mutate(CHR = gsub(pattern = "^chr", replacement = "", do.call(rbind, strsplit(as.character(Locus), split = ':', fixed = TRUE))[,1]),
             POS = as.integer(do.call(rbind, strsplit(as.character(Locus), split = ':', fixed = TRUE))[,2])) %>% 
      subset(R2 >= r2_threshold)
  })
  proxies
}



# rsids_and_chr is a dataframe that contains the following two columns:
# rsid - the rsid of the SNP in question
# CHR - the contig of the SNP in question

# token is an access token to the nci API. 
# if you don't have a token go here: https://ldlink.nci.nih.gov/?tab=apiaccess

get_LD_pairs <- function(rsids_and_chr, population, token) { 
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

# expecting replacement string on the form "A=B,C=D" meaning that
# A will be replaced with B and
# C will be replaced with D

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
