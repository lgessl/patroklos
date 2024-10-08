generate_mock_data <- function(
    n_samples = 10,
    n_genes = 2,
    n_na_in_pheno = 3,
    to_csv = NULL,
    return_type = "data"
){
    # expression matrix
    expr_mat <- matrix(
        sample(1:10, n_samples*n_genes, replace = TRUE),
        nrow = n_samples
    ) |> log()
    rownames(expr_mat) <- stringr::str_c("sample_", 1:n_samples)
    colnames(expr_mat) <- stringr::str_c("gene_", 1:n_genes)
    # to expression tibble
    expr_tbl <- expr_mat |> t() |> tibble::as_tibble(rownames = "gene_id")

    # pheno (n_samples x 5 tibble)
    pheno_tbl <- tibble::tibble(.rows = n_samples)
    pheno_tbl[["patient_id"]] <- rownames(expr_mat)
    pheno_tbl[["progression"]] <- rep(1, n_samples)
    pheno_tbl[["pfs_years"]] <- rep(0, n_samples)
    pheno_tbl[["discrete_var"]] <- sample(1:3, size = n_samples, replace = TRUE)
    pheno_tbl[["ipi_age"]] <- sample(0:1, size = n_samples, replace = TRUE, 
        prob = c(0.27, 0.73))
    pheno_tbl[["abc_gcb"]] <- sample(c("ABC", "GCB", "Unclassified"), 
        size = n_samples, prob = c(0.18, 0.52, 0.30), replace = TRUE)
    pheno_tbl[["continuous_var"]] <- rnorm(n_samples, 10, 10)
    pheno_tbl[["ipi"]] <- sample(1:5, size = n_samples, replace = TRUE)
    pheno_tbl[["split_1"]] <- sample(c("train", "test"), size = n_samples, 
        replace = TRUE)

    # Generate a reasonable time-to-event column as follows:
    # 1. y = x %*% beta + noise (beta, noise ~ N(0, 1))
    # 2. Reduce survival of those samples with IPI = i randomly uniformly 
    # by (11-i)/10 to (10-i)/10 (i = 1, 2, 3, 4, 5)
    # 3. Scale y linearly such that y >= 0 and median is 2
    # 4. Introduce 20% censored samples with censoring time uniformly reduced 
    # by 0-100% of the original survival time 
    # All in all, survival follows a linear model scaled according to the IPI
    x_cont <- cbind(expr_mat, as.matrix(pheno_tbl[, c("continuous_var", "ipi")]))
    colnames(x_cont)[n_genes + 1:2] <- c("continuous_var++", "ipi++")
    x_cat <- dichotomize_tibble(pheno_tbl[, c("discrete_var", "ipi")])
    colnames(x_cat) <- paste0(colnames(x_cat), "++")
    x <- cbind(x_cont, x_cat)
    attr(x, "li_var_suffix") <- "++"
    beta <- rnorm(ncol(x))
    y <- x %*% beta 
    y <- y[, 1]
    for (i in 1:5) {
        y[pheno_tbl[["ipi"]] == i] <- runif(sum(pheno_tbl[["ipi"]] == i), 
            (10-i)/10, (11-i)/10) * y[pheno_tbl[["ipi"]] == i]
    }
    y <- y - quantile(y, 0.76) # 0.76 quantile will become 0
    y <- y/(-min(y)) * 2 # Scale into [-2, \infty)
    y <- y + 2 # 0.76 quantile is 2, minimum is 0
    # y_bin as you can give it to fitter
    y_bin <- as.matrix(as.numeric(y < 2))
    y_bin[1, 1] <- NA
    rownames(y_bin) <- rownames(expr_mat)
    # Censor
    progression <- rep(1, n_samples)
    censored <- sample(seq(n_samples), floor(0.2*n_samples))
    progression[censored] <- 0
    y_uncensored <- y
    y[censored] <- runif(length(censored), 0, 1) * y[censored] 
    y_cox <- cbind(y, progression)
    rownames(y_cox) <- rownames(expr_mat)
    x3y <- list()
    x3y[["x"]] <- x
    x3y[["bin"]] <- binarize_y(y_cox, 1.75, 2)
    x3y[["cox"]] <- censor_y(y_cox, 3)
    x3y[["true"]] <- binarize_y(y_cox, 2, 2)
    if (return_type == "fitter") return(x3y)

    pheno_tbl[["pfs_years"]] <- y
    pheno_tbl[["progression"]] <- progression
    # insert NAs
    pheno_tbl[["ipi"]][1] <- NA
    na_rows <- sample(1:n_samples, n_na_in_pheno, replace = TRUE)
    na_cols <- sample(c("continuous_var", "discrete_var", "ipi"), n_na_in_pheno, 
        replace = TRUE)
    if(n_na_in_pheno > 0){
        for(i in 1:n_na_in_pheno){
            pheno_tbl[na_rows[i], na_cols[i]] <- NA
        }
    }

    if(is.character(to_csv)){
        if(!dir.exists(to_csv)) dir.create(to_csv, recursive = TRUE)
        readr::write_csv(expr_tbl, file.path(to_csv, "expr.csv"))
        readr::write_csv(pheno_tbl, file.path(to_csv, "pheno.csv"))
    }
    if(is.null(to_csv)) to_csv <- "mock_dir"

    data <- Data$new(
        name = "mock",
        directory = to_csv,
        pivot_time_cutoff = 2,
        expr_file = "expr.csv",
        pheno_file = "pheno.csv",
        cohort = "train",
        patient_id_col = "patient_id",
        time_to_event_col = "pfs_years",
        event_col = "progression",
        cohort_col = "split_1",
        benchmark_col = "ipi",
        gene_id_col = "gene_id"
    )
    data$expr_mat <- expr_mat
    data$pheno_tbl <- pheno_tbl
    return(data)
}

apb <- function(
    n_samples,
    fluctuating_availability = TRUE
){
    l <- list(actual = NULL, predicted = NULL, benchmark = NULL)
    for(i in 1:3){
        # Simulate fluctuating availability
        if(fluctuating_availability)
            n_samples <- n_samples + sample(c(-1, 1), size = 1)
        if(i == 1){
            l[[i]] <- sample(c(0, 1), n_samples, replace = TRUE)
        } else {
            l[[i]] <- rnorm(n_samples)
        }
        names(l[[i]]) <- paste0("sample_", 1:n_samples)
        l[[i]][sample(1:n_samples, 1)] <- NA
    }
    return(l)
}

ipi_model <- function(data) {
    ipi_model <- Model$new(
        name = "ipi",
        directory = file.path(data$directory, "ipi"),
        fitter = ptk_zerosum,
        time_cutoffs = 2,
        val_error_fun = neg_prec_with_prev_greater(0.17),
        hyperparams = list(family = "gaussian", alpha = 1, lambda = 0, zeroSum = FALSE, nFold = 2),
        include_from_continuous_pheno = "ipi",
        include_from_discrete_pheno = NULL,
        include_expr = FALSE,
        enable_imputation = FALSE
    )
    data$cohort <- "train"
    ipi_model$fit(data, quiet = TRUE)
    ipi_model$fit_obj$coef[[1]][, 1] <- c(0, 1)
    ipi_model$fit_obj$val_predict[, 1] <- NA
    saveRDS(ipi_model, file.path(dir, "ipi", "model.rds"))
    return(ipi_model)
}