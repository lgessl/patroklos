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
    gl = NULL
){
    if(is.null(gl)) gl <- rep(">", length(col_names))
    stopifnot(all(gl %in% c(">", "<")))
    stopifnot(is.data.frame(tbl))
    stopifnot(is.numeric(cutoffs))
    stopifnot(is.character(gl))
    if(any(c(length(col_names), length(cutoffs)) != length(gl)))
        stop("col_names, cutoffs, and gl must have the same length")
    check_tbl_columns_exist(tbl = tbl, tbl_name = "tbl", col_names = col_names)

    # Extend pheno_tbl
    for(i in seq_along(col_names)){
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
#' the json by hand.
#' @param filename string. The path to the JSON file.
#' @param pheno_tbl tibble. The pheno data. Its format needs to comply with `data_spec` 
#' below.
#' @param expr_tbl tibble. The expression data. Its format needs to comply with `data_spec`
#' below.
#' @param data_spec DataSpec S3 object. Specifications on both `pheno_tbl` and `expr_tbl`.
#' @return The JSON string that was written to the file.
#' @export
write_data_info <- function(
    filename,
    pheno_tbl,
    expr_tbl,
    data_spec
){
    if(file.exists(filename)){
        info_list <- jsonlite::read_json(filename)
    } else {
        info_list <- list(
            "publication" = list(
                "title" = "",
                "author" = "",
                "year" = NA,
                "doi"= "",
                "pubmed id" = NA
            ),
            "data" = list(
                "website" = "",
                "disease type" = "",
                "number of samples" = nrow(pheno_tbl),
                "pheno data" = list(),
                "expression data" = list(),
                "benchmark" = list()
            )
        )
    }

    n_included_in_survival_analysis <- (pheno_tbl[[data_spec$time_to_event_col]] > 0) |> 
        sum(na.rm = TRUE)
    n_high_risk <- ((pheno_tbl[[data_spec$time_to_event_col]] < data_spec$pivot_time_cutoff) & 
        pheno_tbl[[data_spec$event_col]] == 1) |> sum(na.rm = TRUE)
    n_low_risk <- (pheno_tbl[[data_spec$time_to_event_col]] >= 
        data_spec$pivot_time_cutoff) |> sum(na.rm = TRUE)
    prop_high_risk <- n_high_risk / (n_high_risk + n_low_risk)

    pheno_data_list <- list(
        "included in survival analysis" =  n_included_in_survival_analysis,
        "pivot time cutoff" = data_spec$pivot_time_cutoff,
        "number high risk" = n_high_risk,
        "number low risk" = n_low_risk,
        "unknown" = n_included_in_survival_analysis - n_high_risk - n_low_risk,
        "proportion high risk" = prop_high_risk
    )
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
        "name" = data_spec$benchmark_col,
        "reference" = "",
        "precision vs. prevalence" = prec_from_scores(
            pheno_tbl, 
            data_spec
        )
    )
    info_list[["data"]][["pheno data"]] <- pheno_data_list
    info_list[["data"]][["expression data"]] <- expr_data_list
    info_list[["data"]][["benchmark"]] <- benchmark_list

    jsonlite::write_json(info_list, filename, auto_unbox = TRUE, pretty = TRUE,
        dataframe = "columns")
    return(jsonlite::toJSON(info_list, auto_unbox = TRUE, pretty = TRUE, 
        dataframe = "columns"))
}


prec_from_scores <- function(
    pheno_tbl,
    data_spec,
    risk_scores = NULL
){
    if(is.null(risk_scores)){
        if(!is.null(data_spec$benchmark_col)){
            risk_scores <- pheno_tbl[[data_spec$benchmark_col]]
            names(risk_scores) <- pheno_tbl[[data_spec$patient_id_col]]
        } else {
            message("No risk scores provided and no benchmark specified. 
                Returning NULL.")
        }
    }
    model_spec <- list(
        "response_type" = "binary",
        "time_cutoffs" = data_spec$pivot_time_cutoff
    )
    true_risk <- generate_response(
        pheno_tbl = pheno_tbl,
        data_spec = data_spec,
        model_spec = model_spec
    )[, 1]
    tbl <- metric_with_rocr(
        estimate = risk_scores,
        actual = true_risk,
        x_metric = "rpp",
        y_metric = "prec"
    )
    if(is.null(data_spec$benchmark_col)) data_spec$benchmark_col <- "score"
    names(tbl) <- c("rpp", "prec", paste0(data_spec$benchmark_col, " >="))
    tbl <- tbl[, c(3, 1, 2)]
    return(tbl)
}