
#' Find proxy snps to the list provided.
#'
#' function will keep a cache of previously queried proxies in the results dir and only query
#' snps that are not cached. This includes rsids that errored. in that case, an empty file will be
#' stored in the cache directory and the snp will be skipped in a subsequent call. Note that network errors
#' may cause empty files to erroneously indicate that a snp has no rsID.
#'
#' A token for the service is required, please get one at https://ldlink.nci.nih.gov/?tab=apiaccess
#'
#'
#' @name get_proxies
#'
#' @param rsids a list of rsids (of snps) for which proxies are wanted.
#' @param token a token to the ldlink nih service. default will look for an environment variable
#' LDLINK_PROXY. If you need one go here: https://ldlink.nci.nih.gov/?tab=apiaccess
#' @param population specify within which population you want to find proxies. E.g. "CEU", "YRI", etc.
#' @param results_dir A subdirectory (a hash that depends on the population and r2_threshold values) of this directory will be used to
#' cache results
#' @param skip_api A boolean indicating whether to only use the cached results by skipping the API calls.
#' @param r2_threshold The R^2 threshold to use when returning results
#'
#' @export
#'
#'
#' @examples
#'
#' \dontrun{ #because it needs a token.
#' # If one does `Sys.setenv(LDLINK_TOKEN=<Your token here>)` first, this will work
#'
#' get_proxies("rs2495477", pop="CEU", results_dir="derived_data")  #returns one proxy
#' get_proxies("rs373780327", pop="CEU", results_dir="derived_data") #returns no proxies
#'
#' get_proxies(c("rs2495477","rs373780327"), pop="CEU", results_dir="derived_data")  #returns one proxy
#' get_proxies(c(), pop="CEU", results_dir="derived_data")  #returns no proxy (empty list!)
#'
#' # note that thanks to the cache, calling either of these a second time will result in an immediate
#' # return value
#'
#'}
#'
get_proxies <-
  function(rsids,
           token = Sys.getenv("LDLINK_TOKEN"),
           population,
           results_dir = NULL,
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
    # message("listing")
    xfun::in_dir(hashed_subdir_full, {
      raw_rs_files = dir(".", "^rs[0-9]*\\.txt")
    # message("listed")
      
      rsFiles <-
        Filter(function(x)
          file.info(x)$size > 1000, raw_rs_files)
      # message("filtered (size)")
      
      missing_rsids <-
        setdiff(rsids, gsub("\\.txt$", "", raw_rs_files))
      # message(glue::glue("found missing snps (length: {length(missing_rsids)})"))
      
      # this takes time and hits the LDlink API.
      if (!skip_api && length(missing_rsids)) {
        lapply(missing_rsids, function(missing_rsid) {
          print(glue::glue("looking for proxies for {missing_rsid}"))
          myfile <- paste0(missing_rsid, ".txt")
          df_proxy <-
            LDlinkR::LDproxy(
              snp = missing_rsid,
              pop = population,
              r2d = "r2",
              token = token
            )
          if (!(grepl("error", df_proxy[1, 1]))) {
            utils::write.table(
              df_proxy,
              file = myfile,
              append = FALSE,
              quote = FALSE,
              row.names = TRUE,
              sep = "\t"
            )
          } else{
            system(glue::glue("touch {myfile}"))
          }
        })
        
        
        rsFiles <-
          Filter(function(x)
            file.info(x)$size > 1000,
            dir(".", "^rs[0-9]*\\.txt"))
      } else if (length(missing_rsids) > 0) {
        warning(
          glue::glue(
            "Not calling the LDproxy api, but there are some missing rsids in the results dir:{missing_rsids}"
          )
        )
      }
      
      
      
      R2 <-
        Locus <-
        Coord <- Alleles <- RS_Number <- rsid <- query_rsid <- NULL
      # only read snps
      proxies <-Filter(x=rsFiles,f=function(x) {rsid = gsub("\\.txt$", "", x); rsid %in% rsids}) %>%
        plyr::ldply(function(x)
          dplyr::mutate(
            query_rsid = gsub("\\.txt$", "", x),
            utils::read.table(x, sep = "\t")
          ) %>%
            subset(R2 >= 0 &
                     grepl("([ACGT]/[ACGT])", x = Alleles))) %>%
        dplyr::mutate(rsid = RS_Number, Locus = Coord) %>%
        subset(query_rsid %in% rsids) %>%
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
        subset(R2 >= r2_threshold) %>% 
        mutate_cond(condition = rsid == ".", rsid = NA)
    })
    if(nrow(proxies)){
      proxies <-proxies %>% dplyr::select(c("CHR","POS","query_rsid","rsid","Alleles", "R2", "Correlated_Alleles" ,"Locus"))
      assert_proxies(proxies)
    }
    
    proxies
  }


#' Gets the LD between all the pairs of snps (that are in the same contig) from the set provided by
#' the input dataframe
#' 
#' @name get_LD_pairs
#' 
#' @param rsids_and_chr is a dataframe that contains (at least) the following two columns: CHR, rsid
#' @param pop the population codes (or vector thereof) in which the LD is calculated e.g. c("CEU"), c("CEU","YRI")
#' @param token is an access token to the nci API.
#' if you don't have a token go here: https://ldlink.nci.nih.gov/?tab=apiaccess
#'
#'
#' @return a dataframe containing the following columns: CHR, RS_number, variable, value where
#' CHR is the contig on which the snps are
#' RS_number and variable are the two snps being compared and
#' value is the R^2 value between them
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



#' method to remove ambiguous (palindromic and AF too close to 0.5) SNPs and select for the smallest p-value within each "query_rsid":
#' @name choose_best_proxies
#' 
#' @param proxies_and_more the dataframe containing results from get_proxies merged with an outcome GWAS.
#' The Alleles column is assumed to have been already "fixed" according to the "Correlated_Alleles" column (using fix_alleles).
#' @param near_half_threshold The distance from AEF=0.5 that is to be considered "unimbiguous" when a palindromic SNP
#' has it.
#' @param validate a boolean indicating whether to validate the resulting gwas 
#' (TRUE by default) use FALSE for debugging.
#'
#' @return The original dataframe, with palindromic SNPs removed and within each query_rsid only the one with the smallest
#' P value retained. In the case of a tie, the one with the smallest distance to the query SNP (located in POS.exp) will be
#' chosen. In the extremely rare case of a further tie, it will be broken randomly.
#' @export
#'
choose_best_proxies <-
  function(proxies_and_more, near_half_threshold = 0.08, validate=TRUE) {
    proxies_and_more$POS <- proxies_and_more$POS.source
    
    if (validate) {
      assert_probabilities(proxies_and_more)
      assert_gwas(proxies_and_more)
    }
    
    EAF <-
    POS.exp <-
    POS.proxy <-
    Alleles <-
    query_rsid <- 
    maf_near_half <- 
    R2 <- 
    distance_to_query <- NULL
    
    proxies_and_more %>%
      dplyr::mutate(
        maf_near_half = abs(EAF - 0.5) <= near_half_threshold,
        distance_to_query = abs(POS.source - POS.proxy)
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
#' @param data_with_rsids an input dataframe that includes column rsid
#' @param ld_pairs the result from calling get_LD_pairs matrix on the collection of rsids (and their contigs) in data_with_rsids
#' @param r2_threshold the (r^2) threshold to use for considering two snps to be in LD
#'
#' @return a subset of the original dataframe, pruned according to some method....
#' @export
#'
prune_snps <-
  function(data_with_rsids, ld_pairs, r2_threshold = 0.05) {
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
    
    list(rsids = subset(data_with_rsids, rsid %notin% removed),
         removed_rsid = removed)
  }


#' Replaces the alleles in a GWAS according to the replacement string provided.
#'
#' @name replace_alleles
#'
#' @param alleles a list of alleles that will be transformed according to the "rule" 
#' in replacement_string
#' @param replacement_string a set of replacement rules in the format "<CHAR_OLD>=<CHAR_NEW>,..."
#' where the appearances of CHAR_OLD will be replaced by CHAR_NEW. 
#' @param reverse_direction a flag indicating that the replacement should occur in the reverse direction
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
replace_alleles <- function(alleles, replacement_string, reverse_direction=FALSE) {
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
  
  if(reverse_direction){
    replace <- unlist(allele_matrix[1, ])
    pattern <- unlist(allele_matrix[2, ])
  } else {
    replace <- unlist(allele_matrix[2, ])
    pattern <- unlist(allele_matrix[1, ])
  }
    
  ret <- alleles %>% sapply(function(x)
    mgsub::mgsub(
      as.character(x),
      pattern = pattern,
      replacement = replace
    ))
  if (length(unique(ret)) != length(alleles)) {
    warning(glue::glue("Got fewer unique alleles than was provided alleles:{alleles}, replacement_string: {replacement_string}, methinks something went wrong...\n"))
  }
  ret
}

#' Since proxies are provided with a "Correlated_Alleles" column that links the old alleles to the new alleles.
#' This method converts the alleles accordingly
#'
#' @name fix_proxy_alleles
#'
#' @param proxies the dataframe containing (at least) the following columns EA, NEA (alleles), Correlated_Alleles
#' @param reverse_direction a flag indicating that the replacement should occur in the reverse direction
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
fix_proxy_alleles <- function(proxies, reverse_direction = FALSE) {
  . <- NULL
  
  fixed_proxies_with_exposure <- plyr::adply(proxies, 1, function(x) {
    replace_alleles(list(x$EA, x$NEA), 
                    x$Correlated_Alleles,
                    reverse_direction) %>%
    {
      data.frame(EA  = .[1],
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
#' @param prune_r2_threshold a threshold (in LD r^2) for two SNPs to be considered "close"
#'
#' @return a list with the pruned dataframe, the list of removed SNPs, and three plots (if ggplot2 is available)
#' @export
#'
remove_linked_snps <-
  function(input_snps,
           pop,
           token = Sys.getenv("LDLINK_TOKEN"),
           prune_r2_threshold = 0.05) {
    
    variable <- RS_number <- value <- NULL
    ld_pairs_raw <- get_LD_pairs(input_snps, pop, token)
    
    unique_rsids <-
      with(ld_pairs_raw, unique(c(
        as.character(RS_number), as.character(variable)
      )))
    ld_pairs <- ld_pairs_raw %>% dplyr::mutate(
      variable = factor(variable, levels = unique_rsids),
      RS_number = factor(RS_number, levels = unique_rsids)
    )
    
    if (requireNamespace("ggplot2",quietly = TRUE)) {
      p <-
        ggplot2::ggplot(ld_pairs) +  
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + 
          ggplot2::geom_bin2d(stat = 'identity', ggplot2::aes(x = RS_number, y = variable))
      
      raw_ld_plot <-
        p + ggplot2::aes(fill = value) + ggplot2::labs(fill = expression(R ^ 2), parse = TRUE)
      thresholded_ld_plot <-
        p + ggplot2::aes(fill = value > prune_r2_threshold) + 
            ggplot2::labs(fill = expression(R ^ 2 > threshold), parse = TRUE)
    }
    
    list_result <-
      prune_snps(input_snps, ld_pairs, prune_r2_threshold)
    pruned_combined_snps <- list_result$rsids
    removed_rsids <- list_result$removed_rsid
    
    if (requireNamespace("ggplot2",quietly = TRUE)) {
      pruned_ld_plot <- ggplot2::ggplot(ld_pairs) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::geom_bin2d(stat = 'identity', ggplot2::aes(x = RS_number, y = variable)) +
        ggplot2::aes(
          fill = value > prune_r2_threshold,
          alpha = RS_number %in% removed_rsids |
            variable %in% removed_rsids
        ) +
        ggplot2::labs(fill = expression(R ^ 2 > threshold),
             alpha = "variant removed",
             parse = TRUE) +
        ggplot2::scale_alpha_manual(values = c(0.3, 1), breaks = c(TRUE, FALSE))
      
    }
    temp <- list(pruned = pruned_combined_snps,
                 removed_rsids = removed_rsids)
    
    if (requireNamespace("ggplot2",quietly = TRUE)) {
      temp$raw_plot <- raw_ld_plot
      temp$thresholded_plot <- thresholded_ld_plot
      temp$pruned_plot <- pruned_ld_plot
    }
    temp
  }