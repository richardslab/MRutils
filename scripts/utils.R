
library(glue)
library(jsonlite)
library(dplyr)
library(LDlinkR)
library(plyr)
library(reshape2)
library(TwoSampleMR)
library(mgsub)
library(purrr)
library(validate)
library(magrittr)

required_headers <- c("rsid", "CHR", "POS", "P", "beta", "EA", "NEA", "EAF", "SE")


valid_contigs <- c(1:22, "X", "Y")

gwas_rules <- validator(field_format(rsid, "rs*"),
                        CHR %in% valid_contigs,
                        POS > 0,
                        field_format(EA, "[ACGT]", type = "regex"),
                        field_format(NEA, "[ACGT]", type = "regex"),
                        in_range(EAF, 0, 1),
                        is.numeric(beta),
                        is.numeric(SE),
                        in_range(SE, min = 0, Inf,strict = TRUE),
                        in_range(P, 0, 1, strict = TRUE)
)


# it is expected that snps has two columns 'CHR' and 'POS'
# indicating the chromosome and the  (1-based) position of the snp in question.

extract_snps_from_bgzip <-
  function(outcome, snps, chr_col = 1, pos_col = 2, comment_char = "#") {
  #nolint (unused variable)
  region <- snps %>%
    alply(1, function(x) with(x, glue("{CHR}:{POS}-{POS}"))) %>% { 
      do.call(paste, .) 
      }
  
  temp_output <- tempfile(pattern = "subsetted_exposure__",
                          tmpdir = tempdir(),
                          fileext = ".txt")
 
  cmd <- glue("tabix -f -b {pos_col} -c '{comment_char}' -s {chr_col} {outcome} {region} > {temp_output}")
  system(cmd)
  col_names <- read.table(outcome,
                          header = TRUE,
                          nrows = 1,
                          comment.char = "") %>% names()
  col_names[1] <- "CHR"
  outcome_data <- read.table(temp_output, col.names = col_names)
  
  outcome_data
}

chimeric <- c("(A/T)", "(T/A)", "(C/G)", "(G/C)")


# tee <- function(x) {
#   print(x)
#   x
# }


filter_and_write_exposure_data <- function(data,
                                           location_prefix=".",
                                           pvalue_threshold,
                                           rare_threshold) {
  
  data %<>% select(all_of(required_headers))
  
  write.table(data, glue("{location_prefix}exp_extracted_SNPs.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  
  significant_snps <- data %>% subset(P < pvalue_threshold & EAF >= rare_threshold & 1 - EAF >= rare_threshold)
  write.table(significant_snps,
              glue("{location_prefix}exp_significant_SNPs.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  significant_snps
}


#' Title Get rsid name from position
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
      url <- glue(baseURL1)
      Sys.sleep(1)
      read_json(url)$data$spdis[[1]]
    },
    error = function(e) {
      warning("There was an error (1):")
      warning(e)
      warning("Trying to swap ref and alt")
      Sys.sleep(1)
      read_json(glue(baseURL1_swapped))$data$spdis[[1]]
    })
    
    #nolint (unused variable)
    pos <- f$position
    #nolint (unused variable)
    seq_id <- f$seq_id
    
    baseURL2 <- "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/{seq_id}:{pos}:{ref}:{alt}/rsids"
    baseURL2_swapped <- "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/{seq_id}:{pos}:{alt}:{ref}/rsids"
    
    id <- tryCatch({
      url <- glue(baseURL2)
      Sys.sleep(1)
      paste0("rs", read_json(url)$data$rsids[[1]])
    },
    error = function(e){
      warning("There was an error (2):")
      warning(e)
      warning("Trying to swap ref and alt")
      url <- glue(baseURL2_swapped)
      Sys.sleep(1)
      id <- read_json(url)$data$rsids[[1]]
      glue("rs{id}")
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
  unknown_ids <- subset( x = gwas, 
                         subset = is.na(rsid), 
                         select = c(CHR, POS, NEA, EA)) %>% transform(build = build)

  ## this can take time and hits the API multiple times....
  withIds <- adply(unknown_ids, 1, function(x) {
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
  
  gwas_with_ids <- mutate(gwas, rsid = as.character(rsid)) %>%
    merge(subset(rsids, select = c(CHR, POS, rsid)), all.x = T, by = c("CHR","POS")) %>%
    mutate(rsid = adply(do.call(cbind, list(rsid.y, rsid.x)), 1, function(x) first(na.omit(x)))$V1, 
                               rsid.x = NULL, 
                               rsid.y = NULL)
  gwas_with_ids
}


`%notin%` <- Negate(`%in%`)


get_proxies <- function(rsids, token, population, results_dir, skip_api = FALSE, r2_threshold = 0.9) {
  
  xfun::in_dir(results_dir, {
      # this takes time and hits the LDlink API.
    if (!skip_api) {
      LDproxy_batch(rsids, pop = population, r2d = "r2", append = F, token = token)
    }
    
    ## the filter is here because when LDproxy_batch failes, it creates a small file with the error message...
    rsFiles <- Filter(function(x) file.info(x)$size > 1000, 
                   dir(".", "^rs[0-9]*\\.txt"))
    

  # only read snps  
  proxies <- ldply(rsFiles, function(x) mutate(query_rsid = gsub("\\.txt$", "", x), read.table(x, sep = "\t")) %>% 
                    subset(R2 >= 0 & grepl("([ACGT]/[ACGT])", x = Alleles))) %>% 
    mutate(rsid = RS_Number, Locus = Coord) %>% 
    select(c(-"Distance", -"Dprime", -"RegulomeDB", -"Function", -"RS_Number", -"Coord")) %>% 
    mutate(CHR = gsub(pattern = "^chr", replacement = "", do.call(rbind, strsplit(as.character(Locus), split = ':', fixed = TRUE))[,1]),
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
  # call LD matrix per each chromosome
  ddply(rsids_and_chr, .progress = "text", .variables = .(CHR), function(x) {
    if (nrow(x) > 1) {
      LDmatrix(x$rsid, pop = population, token, r2d = "r2") %>%  melt(id.vars = "RS_number")
    } else {
      temp <- data.frame(RS_number = as.character(x$rsid[1]), X = 1)
      names(temp)[2] <- as.character(x$rsid[1])
      temp %>% melt(id.vars = "RS_number")
    }
  })
}

prune_snps <- function(rsids_and_chr, ld_pairs, r2_threshold = 0.05) {
  
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
  
  allele_matrix <- strsplit(as.character(replacement_string), ",", fixed = T)[[1]] %>% 
  {llply(strsplit(x = ., split = "=", fixed = T))} %>% 
    {purrr::transpose(.l = .)} %>% 
    do.call(what = rbind)
  
  replace <- unlist(allele_matrix[2,])
  pattern <- unlist(allele_matrix[1,])
  
  a <- alleles %>% sapply(function(x) mgsub(as.character(x), pattern = pattern, replacement = replace))
  a
}


all(replace_alleles(list("C","G"),"C=A,G=T") == list("A","T"))


find_duplicate_snps <- function(gwas) {
  duplicated_rsids <- gwas[which(duplicated(gwas$rsid)),"rsid"]
  gwas[which(gwas$rsid %in% duplicated_rsids),]
}


get_2smr_results <- function(preprocessed_snps){
  exposure <- read_exposure_data(preprocessed_snps, 
                                 sep = '\t', 
                                 snp_col = "rsid",
                                 beta_col = "beta.exp",
                                 se_col = "SE.exp", 
                                 eaf_col = "EAF.exp", 
                                 effect_allele_col = "EA.exp",
                                 other_allele_col = "NEA.exp",
                                 pval_col = "P.exp")
  
  
  outcome <- read_outcome_data(preprocessed_snps,
                               sep = '\t',
                               snp_col = "rsid",
                               beta_col = "beta.out",
                               se_col = "SE.out",
                               eaf_col = "EAF.out",
                               effect_allele_col = "EA.out",
                               other_allele_col = "NEA.out",
                               pval_col = "P.out")
  
  harmonized = harmonise_data(exposure, outcome)
  regressed <- mr(dat = harmonized)
  
  single_snp_results <- mr_singlesnp(harmonized)
  
  leave_one_out <- mr_leaveoneout(harmonized)
  
  list(data = harmonized,
       regressed = regressed,
       single_snp = single_snp_results,
       heterogeneity = mr_heterogeneity(harmonized),
       pleiotropy = mr_pleiotropy_test(harmonized),
       forest_plot = mr_forest_plot(single_snp_results),
       funnel_plot = mr_funnel_plot(single_snp_results),
       density_plot = mr_density_plot(single_snp_results, regressed),
       scatter_plot = mr_scatter_plot(regressed, harmonized),
       leave_one_out = leave_one_out,
       leave_one_out_plot = mr_leaveoneout_plot(leave_one_out)
  )
}

