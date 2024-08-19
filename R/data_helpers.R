data_initialize <- function(self, private, name, directory, pivot_time_cutoff, 
    expr_file, pheno_file, cohort, patient_id_col, time_to_event_col, event_col, 
    cohort_col, benchmark_col, gene_id_col, imputer){
    stopifnot(is.character(name))
    stopifnot(is.character(directory))
    stopifnot(is.numeric(pivot_time_cutoff))
    stopifnot(is.character(expr_file))
    stopifnot(is.character(pheno_file))
    stopifnot(is.character(cohort))
    stopifnot(is.character(patient_id_col))
    stopifnot(is.character(time_to_event_col))
    stopifnot(is.character(event_col))
    stopifnot(is.character(cohort_col))
    stopifnot(is.character(benchmark_col) || is.null(benchmark_col))
    stopifnot(is.character(gene_id_col))
    stopifnot(is.null(imputer) || is.function(imputer))
    if (!is.null(imputer) && !setequal(names(formals(imputer)), c("x")))
        stop("Imputer must be NULL or a function with a single argument `x`.")

    self$name <- name
    self$directory <- directory
    self$pivot_time_cutoff <- pivot_time_cutoff
    self$expr_file <- expr_file
    self$pheno_file <- pheno_file
    self$cohort <- cohort
    self$patient_id_col <- patient_id_col
    self$time_to_event_col <- time_to_event_col
    self$event_col <- event_col
    self$cohort_col <- cohort_col
    self$benchmark_col <- benchmark_col
    self$gene_id_col <- gene_id_col
    self$imputer <- imputer

    invisible(self)
}

data_read <- function(self, private){
    # extract values from self
    directory <- self$directory
    expr_file <- self$expr_file
    pheno_file <- self$pheno_file
    patient_id_col <- self$patient_id_col
    gene_id_col <- self$gene_id_col

    # read
    self$pheno_tbl <- readr::read_csv(file.path(directory, pheno_file), show_col_types = FALSE)
    expr_tbl <- readr::read_csv(file.path(directory, expr_file), show_col_types = FALSE)

    self$qc_preprocess(expr_tbl = expr_tbl)

    # expression to matrix
    expr_mat <- expr_tbl |>
        dplyr::select(!dplyr::all_of(gene_id_col)) |>
        as.matrix() |> 
        t()
    colnames(expr_mat) <- expr_tbl[[gene_id_col]]
    self$expr_mat <- expr_mat

    invisible(self)
}

data_survival_quantiles <- function(self, private, round_digits) {
    if (is.null(self$pheno_tbl)){
        stop("You need to read in the data first.")
    }
    time_to_event <- self$pheno_tbl[[self$time_to_event_col]]
    progression <- self$pheno_tbl[[self$event_col]]
    order_tte <- order(time_to_event)
    time_to_event <- time_to_event[order_tte]
    progression <- progression[order_tte]
    values  <- unique(time_to_event)

    quantiles <- vapply(values, function(x){
        n_below <- sum(time_to_event < x & progression == 1)
        n_above <- sum(time_to_event >= x)
        round(n_below / (n_below+n_above), digits = round_digits)
        },
        numeric(1)
    )
    quantiles <- round(quantiles, digits = round_digits) 
    n_leading_zeros <- sum(quantiles == "0")
    tbl <- tibble::tibble(quantiles, values)
    names(tbl) <- c("quantile", self$time_to_event_col)
    tbl <- tbl[seq(n_leading_zeros, length(quantiles)), ]
    return(tbl)
}

data_split <- function(self, private, train_prop, save, keep_risk, quiet) {

    stopifnot(train_prop > 0 && train_prop < 1)
    stopifnot(is.logical(save))
    stopifnot(is.logical(keep_risk))
    split_file <- file.path(self$directory, "cohort.rds")
    if (file.exists(split_file)) {
        if (!quiet)
            message("Found a cohort file: ", split_file, ". Loading it and writing ",
                "it into the pheno_tbl attribute.")
        cohort <- readRDS(split_file)
        if (length(cohort) != nrow(self$pheno_tbl) || 
            any(names(cohort) != self$pheno_tbl[[self$patient_id_col]]))
            stop("The patient IDs in the cohort file do not match the patient IDs 
                in the pheno_tbl.")
    } else {
        risk <- rep("na", nrow(self$pheno_tbl))
        if (!is.null(self$pivot_time_cutoff)) {
            risk[self$pheno_tbl[[self$time_to_event_col]] < self$pivot_time_cutoff &
                self$pheno_tbl[[self$event_col]] == 1] <- "high"
            risk[self$pheno_tbl[[self$time_to_event_col]] >= self$pivot_time_cutoff] <- "low"
        }
        risk <- as.factor(risk)
        train_index <- create_data_partition(risk, p = train_prop)
        cohort <- rep("test", nrow(self$pheno_tbl))
        cohort[train_index] <- "train"
        names(cohort) <- self$pheno_tbl[[self$patient_id_col]]
        if (save)
            saveRDS(cohort, split_file)
    }
    self$pheno_tbl[[self$cohort_col]] <- cohort

    invisible(self)
}