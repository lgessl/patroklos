data_initialize <- function(self, private, name, directory, train_prop, 
    pivot_time_cutoff, expr_file, pheno_file, cohort, patient_id_col, 
    time_to_event_col, split_col_prefix, event_col, benchmark_col, gene_id_col){
    
    if(!is.null(cohort))
        cohort <- match.arg(cohort, c("train", "test"))
    stopifnot(is.character(name))
    stopifnot(is.character(directory))
    stopifnot(is.numeric(train_prop))
    stopifnot(is.numeric(pivot_time_cutoff))
    stopifnot(is.character(expr_file))
    stopifnot(is.character(pheno_file))
    stopifnot(is.null(cohort) || is.character(cohort))
    stopifnot(is.character(patient_id_col))
    stopifnot(is.character(time_to_event_col))
    stopifnot(is.character(event_col))
    stopifnot(is.character(benchmark_col) || is.null(benchmark_col))
    stopifnot(is.character(gene_id_col))
    stopifnot(is.character(split_col_prefix))

    self$name <- name
    self$directory <- directory
    self$train_prop <- train_prop
    self$pivot_time_cutoff <- pivot_time_cutoff
    self$expr_file <- expr_file
    self$pheno_file <- pheno_file
    self$cohort <- cohort
    self$patient_id_col <- patient_id_col
    self$time_to_event_col <- time_to_event_col
    self$event_col <- event_col
    self$benchmark_col <- benchmark_col
    self$gene_id_col <- gene_id_col
    self$split_col_prefix <- split_col_prefix

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
    files <- c(expr_file, pheno_file)
    tbls <- list()
    for (i in 1:length(files)){
        full_path <- file.path(directory, files[i])
        tbls[[i]] <- readr::read_csv(full_path, show_col_types = FALSE)
    }
    expr_tbl <- tbls[[1]]
    pheno_tbl <- tbls[[2]]

    # check if identifier columns are there
    if(is.null(expr_tbl[[gene_id_col]])){
        stop("There is no column named ", gene_id_col, " in ", expr_file)
    }
    if(is.null(pheno_tbl[[patient_id_col]])){
        stop("There is no column named ", patient_id_col, " in ", pheno_file)
    }

    # expression to matrix
    gene_names <- expr_tbl[[gene_id_col]]
    if(!elements_unique(gene_names)){
        stop("Column ", gene_id_col, " in ", expr_file, " holds duplicate entries.")
    }
    expr_mat <- expr_tbl |>
        dplyr::select(!dplyr::all_of(gene_id_col)) |>
        as.matrix() |> 
        t()
    colnames(expr_mat) <- gene_names

    # pheno: move patient ids into first column
    patient_ids <- pheno_tbl[[patient_id_col]]
    if(!elements_unique(patient_ids)){
        stop("Column ", patient_id_col, " in ", pheno_file, " holds duplicate entries.")
    }
    pheno_tbl <- pheno_tbl |>
        dplyr::relocate(dplyr::all_of(patient_id_col))

    self$expr_mat <- expr_mat
    self$pheno_tbl <- pheno_tbl

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