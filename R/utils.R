


filter_and_write_exposure_data <- function(data,
                                           location_prefix=".",
                                           pvalue_threshold,
                                           rare_threshold) {
  
  P <- EAF <- NULL
  
  data <- data %>% dplyr::select(dplyr::all_of(required_headers))
  
  utils::write.table(data, glue::glue("{location_prefix}exp_extracted_SNPs.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  
  significant_snps <- data %>% subset(P < pvalue_threshold & EAF >= rare_threshold & 1 - EAF >= rare_threshold)
  utils::write.table(significant_snps,
              glue::glue("{location_prefix}exp_significant_SNPs.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  significant_snps
}



find_duplicate_snps <- function(gwas) {
  duplicated_rsids <- gwas[which(duplicated(gwas$rsid)),"rsid"]
  gwas[which(gwas$rsid %in% duplicated_rsids),]
}


get_2smr_results <- function(preprocessed_snps){
  exposure <- TwoSampleMR::read_exposure_data(preprocessed_snps, 
                                 sep = '\t', 
                                 snp_col = "rsid",
                                 beta_col = "beta.exp",
                                 se_col = "SE.exp", 
                                 eaf_col = "EAF.exp", 
                                 effect_allele_col = "EA.exp",
                                 other_allele_col = "NEA.exp",
                                 pval_col = "P.exp")
  
  
  outcome <- TwoSampleMR::read_outcome_data(preprocessed_snps,
                               sep = '\t',
                               snp_col = "rsid",
                               beta_col = "beta.out",
                               se_col = "SE.out",
                               eaf_col = "EAF.out",
                               effect_allele_col = "EA.out",
                               other_allele_col = "NEA.out",
                               pval_col = "P.out")
  
  harmonized = TwoSampleMR::harmonise_data(exposure, outcome)
  regressed <- TwoSampleMR::mr(dat = harmonized)
  
  single_snp_results <- TwoSampleMR::mr_singlesnp(harmonized)
  
  leave_one_out <- TwoSampleMR::mr_leaveoneout(harmonized)
  
  list(data = harmonized,
       regressed = regressed,
       single_snp = single_snp_results,
       heterogeneity = TwoSampleMR::mr_heterogeneity(harmonized),
       pleiotropy = TwoSampleMR::mr_pleiotropy_test(harmonized),
       forest_plot = TwoSampleMR::mr_forest_plot(single_snp_results),
       funnel_plot = TwoSampleMR::mr_funnel_plot(single_snp_results),
       density_plot = TwoSampleMR::mr_density_plot(single_snp_results, regressed),
       scatter_plot = TwoSampleMR::mr_scatter_plot(regressed, harmonized),
       leave_one_out = leave_one_out,
       leave_one_out_plot = TwoSampleMR::mr_leaveoneout_plot(leave_one_out)
  )
}

