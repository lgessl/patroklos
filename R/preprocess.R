# Helpful functions for preprocessing

#' @title Discretize columns of a tibble and append the resulting columns
#' @description Discretize the columns in a pheno tibble that belong to the features of
#' the International Prognostic Index (IPI) into binary indicator columns by comparing them
#' to cutoff values and append them to the pheno tibble.
#' @param tbl tibble.
#' @param col_names character vector. Discretize the columns in `tbl` that have these
#' names.
#' @param cutoffs numeric vector of the same length as character vector. The cutoff values
#' corresponding to the columns in `col_names`.
#' @param gl character vector of the same length as `col_names` holding ">" and "<"
#' entries. Discretize via greater (`> cutoffs`) or less (`< cutoffs`). Default is `NULL`,
#' in which case all columns are discretized via `> cutoffs`.
#' @return `tbl` with discretized columns as numeric columns appended.
#' @export
discretize_tbl_cols <- function(
    tbl,
    col_names,
    cutoffs,
    gl = NULL) {
    if (is.null(gl)) gl <- rep(">", length(col_names))
    stopifnot(all(gl %in% c(">", "<")))
    stopifnot(is.data.frame(tbl))
    stopifnot(is.numeric(cutoffs))
    stopifnot(is.character(gl))
    if (any(c(length(col_names), length(cutoffs)) != length(gl))) {
        stop("col_names, cutoffs, and gl must have the same length")
    }
    check_tbl_columns_exist(tbl = tbl, tbl_name = "tbl", col_names = col_names)

    # Extend pheno_tbl
    for (i in seq_along(col_names)) {
        old_col <- tbl[[col_names[i]]]
        new_col <- ifelse(
            rep(gl[i] == ">", length(old_col)),
            tbl[[col_names[i]]] > cutoffs[i],
            tbl[[col_names[i]]] < cutoffs[i]
        )
        tbl[[paste0(col_names[i], gl[i], cutoffs[i])]] <- as.numeric(new_col)
    }
    return(tbl)
}


#' @title Write data info to a JSON file
#' @description Automatically generate info about the data and write it to a JSON file.
#' If the file already exists, the automatically inferable info into the file. If the
#' file is not present yet, create a sceleton, write the the info one can automtically
#' infer from the data into the json. You can then fill the empty fields in by writing
#' the json by hand. Manullay alrady filled out info will not be overwritten.
#' @param filename string. The path to the JSON file.
#' @param data Data S3 object. Specifications on both `pheno_tbl` and `expr_tbl`.
#' @param expr_tbl tibble. The expression data. Its format needs to comply with `data`
#' below.
#' @return The JSON string that was written to the file.
#' @export
write_data_info <- function(
    filename,
    data,
    expr_tbl
    ) {
    found_file <- FALSE
    if (file.exists(filename)) {
        info_list <- jsonlite::read_json(filename)
        found_file <- TRUE
    } else {
        info_list <- list(
            "publication" = list(
                "title" = "",
                "author" = "",
                "year" = NA,
                "doi" = "",
                "pubmed id" = NA
            ),
            "data" = list(
                "website" = "",
                "disease type" = "",
                "number of samples" = nrow(data$pheno_tbl),
                "pheno data" = list(),
                "expression data" = list(),
                "benchmark" = list()
            )
        )
    }

    risk_stat_list <- list()
    for (cohort in c(".", unique(data$pheno_tbl[[data$cohort_col]]))) {
        pheno_tbl <- data$pheno_tbl[stringr::str_detect(
            data$pheno_tbl[[data$cohort_col]],
            cohort
        ), ]
        n_included_in_survival_analysis <- (pheno_tbl[[data$time_to_event_col]] > 0) |>
            sum(na.rm = TRUE)
        n_high_risk <- ((pheno_tbl[[data$time_to_event_col]] < data$pivot_time_cutoff) &
            pheno_tbl[[data$event_col]] == 1) |> sum(na.rm = TRUE)
        n_low_risk <- (pheno_tbl[[data$time_to_event_col]] >=
            data$pivot_time_cutoff) |> sum(na.rm = TRUE)
        prop_high_risk <- n_high_risk / (n_high_risk + n_low_risk)
        risk_stat_list[[cohort]] <- list(
            "included in survival analysis" = n_included_in_survival_analysis,
            "pivot time cutoff" = data$pivot_time_cutoff,
            "number high risk" = n_high_risk,
            "number low risk" = n_low_risk,
            "unknown" = n_included_in_survival_analysis - n_high_risk - n_low_risk,
            "proportion high risk" = prop_high_risk
        )
    }
    names(risk_stat_list)[1] <- "all"

    expr_data_list <- list(
        "technology" = "",
        "platform" = "",
        "read length" = "",
        "number of genes" = nrow(expr_tbl),
        "gene symbols" = "",
        "primary processing" = "",
        "normalization" = ""
    )
    benchmark_list <- list(
        "name" = data$benchmark_col,
        "reference" = "",
        "performance" = list()
    )
    for (cohort in c(".", unique(data$pheno_tbl[[data$cohort_col]]))) {
        pheno_tbl <- data$pheno_tbl[stringr::str_detect(
            data$pheno_tbl[[data$cohort_col]],
            cohort
        ), ]
        ch_data <- data$clone()
        ch_data$pheno_tbl <- pheno_tbl
        benchmark_list[["performance"]][[cohort]] <- prec_from_scores(ch_data)
    }
    names(benchmark_list[["performance"]])[1] <- "all"
    info_list[["data"]][["risk stats"]] <- risk_stat_list
    if (!found_file) {
        info_list[["data"]][["expression data"]] <- expr_data_list
        info_list[["data"]][["benchmark"]] <- benchmark_list
    } else {
        info_list[["data"]][["expression data"]][["number of genes"]] <-
            expr_data_list[["number of genes"]]
        info_list[["data"]][["benchmark"]][-2] <- benchmark_list[-2]
    }

    jsonlite::write_json(info_list, filename,
        auto_unbox = TRUE, pretty = TRUE,
        dataframe = "columns"
    )
    invisible(jsonlite::toJSON(info_list,
        auto_unbox = TRUE, pretty = TRUE,
        dataframe = "columns"
    ))
}


prec_from_scores <- function(
    data,
    risk_scores = NULL) {
    if (is.null(data$pheno_tbl)) {
        stop("Data must have `pheno_tbl` read in")
    }
    if (is.null(risk_scores)) {
        if (!is.null(data$benchmark_col)) {
            risk_scores <- data$pheno_tbl[[data$benchmark_col]]
            names(risk_scores) <- data$pheno_tbl[[data$patient_id_col]]
        } else {
            message("No risk scores provided and no benchmark specified.
                Returning NULL.")
        }
    }
    model <- list(
        "time_cutoffs" = data$pivot_time_cutoff
    )
    y_cox <- as.matrix(data$pheno_tbl[, c(data$time_to_event_col, data$event_col)])
    rownames(y_cox) <- data$pheno_tbl[[data$patient_id_col]]
    true_risk <- binarize_y(y_cox,
        time_cutoff = data$pivot_time_cutoff,
        pivot_time_cutoff = data$pivot_time_cutoff
    )[, 1]
    ea <- intersect_by_names(risk_scores, true_risk, rm_na = c(TRUE, TRUE))
    tbl <- metric_with_rocr(
        estimate = ea[[1]],
        actual = ea[[2]],
        x_metric = "rpp",
        y_metric = "prec"
    )
    if (is.null(data$benchmark_col)) data$benchmark_col <- "score"
    names(tbl) <- c("rpp", "prec", paste0(data$benchmark_col, " >="))
    tbl <- tbl[, c(3, 1, 2)]
    return(tbl)
}
