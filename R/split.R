#' @title Split a data set into train and test data multiple times
#' @description Split a data set into train and test data sets by noting down 
#' the affiliation to the train or test cohort in the pheno data for every split
#' Optionally, you can preserve the ratio of high- and low-risk patients. Ensure 
#' that every model in a list of `Model`s gets enough splits.
#' @param pheno_tbl A tibble holding the pheno data (see `Data()`
#' for details).
#' @param data A `Data` object referring to `expr_tbl` and `pheno_tbl`.
#' @param model_list A list of `Model` objects. Ensure that every model
#' gets the splits it needs.
#' @return `pheno_tbl` with all required splits.
ensure_splits <- function(
    pheno_tbl,
    data,
    model_list
){
    # Find split indices to generate new splits for
    wanted_splits <- sapply(
        model_list, 
        function(model) model$split_index
    ) |> unlist() |> as.integer()
    present_splits <- stringr::str_extract(
        names(pheno_tbl),
        stringr::str_c(data$split_col_prefix, "[0-9]+")
    ) |> stringr::str_extract("[0-9]+$") |> as.integer()
    new_splits <- setdiff(wanted_splits, present_splits)

    # Extract from data
    time_col <- data$time_to_event_col
    event_col <- data$event_col
    pivot_time_cutoff <- data$pivot_time_cutoff
    train_prop <- data$train_prop

    # Generate new splits
    for(i in seq_along(new_splits)){
        if(!is.null(data$pivot_time_cutoff)){
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
        pheno_tbl[[stringr::str_c(data$split_col_prefix, new_splits[i])]] <- 
            affiliation
    }

    # Write back to pheno csv
    if(length(new_splits) > 0){
        pheno_file <- file.path(data$directory, data$pheno_file)
        readr::write_csv(pheno_tbl, pheno_file)
    }
    
    return(pheno_tbl)
}