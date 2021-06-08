







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
#'
get_proxies <-
  function(rsids,
           token = Sys.getenv("LDLINK_TOKEN"),
           population,
           results_dir,
           skip_api = FALSE,
           r2_threshold = 0.9) {
    config = list(pop = population, r2_thresh = r2_threshold)
    hashed_subdir = digest::digest(config)
    hashed_subdir_full = glue::glue("{results_dir}/{hashed_subdir}")
    
    if (!dir.exists(hashed_subdir_full)) {
      dir.create(hashed_subdir_full)
      xfun::in_dir(hashed_subdir_full, {
        yaml::write_yaml(config, "config.yaml")
      })
    }
    
    xfun::in_dir(hashed_subdir_full, {
      raw_rs_files = dir(".", "^rs[0-9]*\\.txt")
      
      rsFiles <-
        Filter(function(x)
          file.info(x)$size > 1000, raw_rs_files)
      
      missing_rsids <-
        setdiff(rsids, gsub("\\.txt$", "", raw_rs_files))
      
      # this takes time and hits the LDlink API.
      if (length(missing_rsids) > 0) {
        if (!skip_api) {
          LDlinkR::LDproxy_batch(
            missing_rsids,
            pop = population,
            r2d = "r2",
            append = F,
            token = token
          )
          rsFiles <-
            Filter(function(x)
              file.info(x)$size > 1000,
              dir(".", "^rs[0-9]*\\.txt"))
        } else {
          warning(
            glue::glue(
              "Not calling the LDproxy api, but there are some missing rsids in the results dir:{missing_rsids}"
            )
          )
        }
      }
      
      
      R2 <-
        Locus <-
        Coord <- Alleles <- RS_Number <- rsid <- query_rsid <- NULL
      # only read snps
      proxies <-
        plyr::ldply(rsFiles, function(x)
          dplyr::mutate(
            query_rsid = gsub("\\.txt$", "", x),
            utils::read.table(x, sep = "\t")
          ) %>%
            subset(R2 >= 0 &
                     grepl("([ACGT]/[ACGT])", x = Alleles))) %>%
        dplyr::mutate(rsid = RS_Number, Locus = Coord) %>%
        subset(query_rsid %in% rsids) %>%
        dplyr::select(c(
          -"Distance",
          -"Dprime",
          -"RegulomeDB",
          -"Function",
          -"RS_Number",
          -"Coord",
          -"MAF"
        )) %>%
        dplyr::mutate(
          CHR = gsub(
            pattern = "^chr",
            replacement = "",
            do.call(rbind, strsplit(
              as.character(Locus), split = ':', fixed = TRUE
            ))[, 1]
          ),
          POS = as.integer(do.call(
            rbind, strsplit(as.character(Locus), split = ':', fixed = TRUE)
          )[, 2])
        ) %>%
        subset(R2 >= r2_threshold) %>% mutate_cond(condition = rsid == ".", rsid =
                                                     NA)
    })
    assert_proxies(proxies)
    proxies
  }



#' Gets the LD between all the pairs of snps (that are in the same contig) from the set provbided by
#' the input dataframe
#'
#' @param rsids_and_chr is a dataframe that contains (at least) the following two columns: CHR, rsid
#' @param pop the population codes (or vector thereof) in which the LD is calculated e.g. c("CEU"), c("CEU","YRI")
#' @param token is an access token to the nci API.
#' if you don't have a token go here: https://ldlink.nci.nih.gov/?tab=apiaccess
#'
#'
#' @return
#' 
#' @export
#'
#' @examples
#'
#' \dontrun{
#'    # returns DF with 4 rows
#'    get_LD_pairs(data.frame(CHR="chr9", 
#'                            rsid=c("rs10760259","rs635634")),
#'                            pop=c("CEU"), 
#'                            token=Sys.getenv("LDLINK_TOKEN"))
#' }
#'
get_LD_pairs <- function(rsids_and_chr, pop, token) {
  assert_rsids(rsids_and_chr)
  assert_chr(rsids_and_chr)
  assertthat::assert_that(
    nchar(token) > 1,
    msg = glue::glue(
      "Token provided ({token}), is invalid, please obtain one at https://ldlink.nci.nih.gov/?tab=apiaccess"
    )
  )
  CHR <- NULL
  # call LD matrix per each chromosome
  plyr::ddply(rsids_and_chr, .progress = if(interactive()){"text"} else {"none"}, .variables = "CHR", function(x) {
    if (nrow(x) > 1) {
      suppressMessages(LDlinkR::LDmatrix(x$rsid, pop = pop, token, r2d = "r2") %>%
        reshape2::melt(id.vars = "RS_number"))
    } else {
      temp <- data.frame(RS_number = as.character(x$rsid[1]), X = 1)
      names(temp)[2] <- as.character(x$rsid[1])
      temp %>% reshape2::melt(id.vars = "RS_number")
    }
  })
}



#' method to remove ambiguous (palindromic and AF too close to 0.5) snps and select for the smallest pvalue within each "query_rsid":
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
choose_best_proxies <-
  function(proxies_and_more, near_half_threshold = 0.08) {
    proxies_and_more$POS <- proxies_and_more$POS.exp
    assert_probabilities(proxies_and_more)
    assert_gwas(proxies_and_more)
    
    EAF <-
      POS.exp <-
      POS.proxy <-
      Alleles <-
      query_rsid <- maf_near_half <- R2 <- distance_to_query <- NULL
    
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
      dplyr::mutate(
        maf_near_half = NULL,
        distance_to_query = NULL,
        POS = NULL
      )
  }


#' Prunes a collection of SNPs provided in a dataframe.
#' The method is as follows:
#'
#' @name prune_snps
#'
#' @param data_with_rsids an input dataframe that includes columns rsid
#' @param ld_pairs the result from calling get_LD_pairs matrix on the collection of rsids (and their contigs) in rsids_and_chr
#' @param r2_threshold the (r^2) threshold to use for considering two snps to be in LD
#'
#' @return a subset of the original dataframe, pruned according to some method....
#' @export
#'
prune_snps <-
  function(rsids_and_chr, ld_pairs, r2_threshold = 0.05) {
    RS_number <- rsid <- value <- variable <- NULL
    
    LDpairs_culled <- ld_pairs
    pairs <-
      subset(
        LDpairs_culled,
        value >= r2_threshold &
          as.character(RS_number) != as.character(variable)
      )
    removed <- c()
    while (nrow(pairs) > 0) {
      ##TODO: make sure to keep snps with smaller P value.....
      to_remove <-
        names(which.max(table(c(
          as.character(pairs$RS_number),
          as.character(pairs$variable)
        ))))
      # to_remove <- pairs[which.max(pairs$value), "variable"] %>% as.character()
      removed <- c(removed, to_remove)
      
      LDpairs_culled <-
        subset(LDpairs_culled,
               variable != to_remove & RS_number != to_remove)
      pairs <-
        subset(
          LDpairs_culled,
          value >= r2_threshold &
            as.character(RS_number) != as.character(variable)
        )
    }
    
    list(rsids = subset(rsids_and_chr, rsid %notin% removed),
         removed_rsid = removed)
  }


#' Replaces the alleles in a gwas according to the replacement string provided.
#'
#' @name replace_alleles
#'
#' @param alleles
#' @param replacement_String
#'
#' @keywords internal
#' @export
#' @examples
#' replace_alleles(list("C", "G"), "C=A,G=T") # c("A", "T")
#' replace_alleles(list("C", "G"), "C=C,G=T") # c("C", "T")
#' replace_alleles(list("C", "G"), "C=A,G=G") # c("A", "G")
#' replace_alleles(list("C", "G"), "C=T,G=A") # c("T", "A")
#' replace_alleles(list("C", "T"), "C=A,T=T") # c("A", "T")
#' replace_alleles(list("C", "T"), "C=A,T=C") # c("A", "C")
#' replace_alleles(list("C", "T"), "C=T,T=C") # c("T", "C")
#'
replace_alleles <- function(alleles, replacement_string) {
  . <- NULL
  allele_matrix <-
    strsplit(as.character(replacement_string), ",", fixed = T)[[1]] %>%
    {
      plyr::llply(strsplit(
        x = .,
        split = "=",
        fixed = T
      ))
    } %>%
    {
      purrr::transpose(.l = .)
    } %>%
    do.call(what = rbind)
  
  replace <- unlist(allele_matrix[2, ])
  pattern <- unlist(allele_matrix[1, ])
  
  ret <- alleles %>% sapply(function(x)
    mgsub::mgsub(
      as.character(x),
      pattern = pattern,
      replacement = replace
    ))
  if (length(unique(ret)) != length(alleles)) {
    warning("Got fewer unique alleles than was provided, methinks something went wrong...")
  }
  ret
}

#' Since proxies are provided with a "Correlated_Alleles" column that links the old alleles to the new alleles.
#' This method converts the alleles accordingly
#'
#' @name fix_proxy_alleles
#'
#' @param proxies the dataframe containing (at least) the following columns EA, NEA (alleles), Correlated_Alleles
#'
#' @return the same dataframe with the alleles in EA and NEA replaced according to Correlated_Alleles column
#'
#' @export
#'
#' @examples
#'
#' demo_data
#' fix_proxy_alleles(data.frame(demo_data,
#'                              Correlated_Alleles=c("G=G,A=T",
#'                                                   "C=G,G=C",
#'                                                   "A=G,G=T",
#'                                                   "C=T,T=A",
#'                                                   "C=C,A=T",
#'                                                   "A=T,G=C")))
#'
fix_proxy_alleles <- function(proxies) {
  fixed_proxies_with_exposure <- plyr::adply(proxies, 1, function(x) {
    replace_alleles(list(x$EA, x$NEA), x$Correlated_Alleles) %>%
    {
      # print(.[1]);
      # print(.[2]);
      data.frame(EA = .[1],
                 NEA = .[2])
    }
  })
  
  fixed_proxies_with_exposure
}


#' Method to select snps so that they are not in LD.
#'
#' @param input_snps a dataframe
#' @param pop the population codes (or vector thereof) in which the LD is calculated e.g. c("CEU"), c("CEU","YRI")
#' @param token is an access token to the nci API.
#' if you don't have a token go here: https://ldlink.nci.nih.gov/?tab=apiaccess
#' @param prune_r2_threshold a threshold (in LD r^2) for two snps to be considered "close"
#'
#' @return a list with the pruned dataframe, the list of removed snps, and three plots (if ggplot2 is available)
#' @export
#'
remove_linked_snps <-
  function(merged_and_combined,
           pop,
           token = Sys.getenv("LDLINK_TOKEN"),
           prune_r2_threshold = 0.05) {
    ld_pairs_raw <- get_LD_pairs(merged_and_combined, pop, token)
    
    unique_rsids <-
      with(ld_pairs_raw, unique(c(
        as.character(RS_number), as.character(variable)
      )))
    ld_pairs <- ld_pairs_raw %>% dplyr::mutate(
      variable = factor(variable, levels = unique_rsids),
      RS_number = factor(RS_number, levels = unique_rsids)
    )
    
    if (require(ggplot2)) {
      p <-
        ggplot(ld_pairs) +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_bin2d(stat = 'identity', aes(x = RS_number, y = variable))
      
      raw_ld_plot <-
        p + aes(fill = value) + labs(fill = expression(R ^ 2), parse = TRUE)
      thresholded_ld_plot <-
        p + aes(fill = value > params$prune_r2_threshold) + labs(fill = expression(R ^
                                                                                     2 > threshold), parse = TRUE)
    }
    
    list_result <-
      prune_snps(merged_and_combined, ld_pairs, prune_r2_threshold)
    pruned_combined_snps <- list_result$rsids
    removed_rsids <- list_result$removed_rsid
    
    if (require(ggplot2)) {
      pruned_ld_plot <- ggplot(ld_pairs) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        geom_bin2d(stat = 'identity', aes(x = RS_number, y = variable)) +
        aes(
          fill = value > params$prune_r2_threshold,
          alpha = RS_number %in% removed_rsids |
            variable %in% removed_rsids
        ) +
        labs(fill = expression(R ^ 2 > threshold),
             alpha = "variant removed",
             parse = TRUE) +
        scale_alpha_manual(values = c(0.3, 1), breaks = c(TRUE, FALSE))
      
    }
    temp <- list(pruned = pruned_combined_snps,
                 removed_rsids = removed_rsids)
    
    if (require(ggplot2)) {
      temp$raw_plot <- raw_ld_plot
      temp$thresholded_plot <- thresholded_ld_plot
      temp$pruned_plot <- pruned_ld_plot
    }
    temp
  }