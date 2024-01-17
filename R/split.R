#' @title Split a data set into train and test data multiple times
#' @description Split a data set into train and test data sets by noting down 
#' the affiliation to the train or test cohort in the pheno data for every split
#' Optionally, you can preserve the ratio of high- and low-risk patients. Ensure 
#' that every model in a list of `ModelSpec`s gets enough splits.
#' @param pheno_tbl A tibble holding the pheno data (see `DataSpec()`
#' for details).
#' @param data_spec A `DataSpec` object referring to `expr_tbl` and `pheno_tbl`.
#' @param model_spec_list A list of `ModelSpec` objects. Ensure that every model
#' gets the splits it needs.
#' @return `pheno_tbl` with all required splits.
ensure_splits <- function(
    pheno_tbl,
    data_spec,
    model_spec_list
){
    # Find split indices to generate new splits for
    wanted_splits <- sapply(
        model_spec_list, 
        function(model_spec) model_spec$split_index
    ) |> unlist() |> as.integer()
    present_splits <- stringr::str_extract(
        names(pheno_tbl),
        stringr::str_c(data_spec$split_col_prefix, "[0-9]+")
    ) |> stringr::str_extract("[0-9]+$") |> as.integer()
    new_splits <- setdiff(wanted_splits, present_splits)

    # Extract from data_spec
    time_col <- data_spec$time_to_event_col
    event_col <- data_spec$event_col
    pivot_time_cutoff <- data_spec$pivot_time_cutoff
    train_prop <- data_spec$train_prop

    # Generate new splits
    for(i in seq_along(new_splits)){
        if(!is.null(data_spec$pivot_time_cutoff)){
            risk <- rep("na", nrow(pheno_tbl))
            risk[pheno_tbl[[time_col]] < pivot_time_cutoff & 
                pheno_tbl[[event_col]] == 1] <- "high"
            risk[pheno_tbl[[time_col]] >= pivot_time_cutoff] <- "low"
            risk <- as.factor(risk)
            train_index <- create_data_partition(risk, p = train_prop)
        } else {
            n_train <- round(nrow(pheno_tbl) * train_prop)
            train_index <- sample(1:nrow(pheno_tbl), n_train, replace = FALSE)
        }
        affiliation <- rep("test", nrow(pheno_tbl))
        affiliation[train_index] <- "train"
        pheno_tbl[[stringr::str_c(data_spec$split_col_prefix, new_splits[i])]] <- 
            affiliation
    }

    # Write back to pheno csv
    if(length(new_splits) > 0){
        pheno_fname <- file.path(data_spec$directory, data_spec$pheno_fname)
        readr::write_csv(pheno_tbl, pheno_fname)
    }
    
    return(pheno_tbl)
}