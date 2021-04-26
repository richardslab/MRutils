
# it is expected that snps has two columns 'CHR' and 'POS'
# indicating the chromosome and the  (1-based) position of the snp in question.

extract_snps_from_bgzip <-
  function(outcome, snps, chr_col = 1, pos_col = 2, comment_char = "#") {
    #nolint (unused variable)
    region <- snps %>%
      plyr::alply(1, function(x) with(x, glue("{CHR}:{POS}-{POS}"))) %>% { 
        do.call(paste, plyr::.) 
      }
    
    temp_output <- tempfile(pattern = "subsetted_exposure__",
                            tmpdir = tempdir(),
                            fileext = ".txt")
    
    cmd <- glue::glue("tabix -f -b {pos_col} -c '{comment_char}' -s {chr_col} {outcome} {region} > {temp_output}")
    system(cmd)
    col_names <- utils::read.table(outcome,
                            header = TRUE,
                            nrows = 1,
                            comment.char = "") %>% names()
    col_names[1] <- "CHR"
    outcome_data <- utils::read.table(temp_output, col.names = col_names)
    
    outcome_data
  }
