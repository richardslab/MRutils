
required_headers <- c("rsid", "CHR", "POS", "P", "beta", "EA", "NEA", "EAF", "SE")

valid_contigs <- c(1:22, "X", "Y")

gwas_rules <- validate::validator(field_format(rsid, "rs*"),
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


palindromic <- c("(A/T)", "(T/A)", "(C/G)", "(G/C)")

