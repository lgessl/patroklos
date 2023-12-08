#' @title Split a data set into train and test data
#' @description Split a data set into train and test data sets, (optionally) 
#' preserving the ratio of high- and low-risk patients.
#' @param expr_tbl A tibble holding the expression data (see `DataSpec()`
#' for details).
#' @param pheno_tbl A tibble holding the pheno data (see `DataSpec()`
#' for details).
#' @param data_spec A `DataSpec` object referring to `expr_tbl` and `pheno_tbl`.
#' @param train_prop numeric. The proportion of samples to be included in the
#' train data set.
#' @param pfs_cut numeric in [0, 1] or NULL. The cut-off for progression-free survival 
#' (PFS) to distinguish between high- and low-risk patients. Default is `NULL`, but must 
#' be specified if `based_on_pfs_cut` is `TRUE`.
#' @param based_on_pfs_cut logical. If `TRUE`, split the data set into high-risk, low-risk,
#' NA (where PFS < `pfs_cut` cannot be inferred from `pheno_tbl`) patients based on `pfs_cut`
#' and preserve the ratio of these three groups in the resulting train and test data sets.
#' Default is `FALSE`, in which case train and test samples are drawn randomly without 
#' any further constraints.
#' @param quiet logical. If `TRUE`, suppress all messages. Default is `FALSE`. 
#' @export
split_dataset <- function(
    expr_tbl,
    pheno_tbl,
    data_spec,
    train_prop,
    pfs_cut = NULL,
    based_on_pfs_cut = FALSE,
    quiet = FALSE
){
    check_consistent_patient_ids(
        stage = "preprocessing", 
        expr_tbl, 
        pheno_tbl, 
        data_spec
    )
    pfs_col <- data_spec$pfs_col
    progression_col <- data_spec$progression_col

    # Split indices
    if(based_on_pfs_cut){
        if(is.null(pfs_cut)){
            stop("pfs_cut must be specified if based_on_pfs_cut is TRUE")
        }
        risk_group <- rep("na", nrow(pheno_tbl))
        risk_group[pheno_tbl[[pfs_col]] < pfs_cut & pheno_tbl[[progression_col]] == 1] <- "high"
        risk_group[pheno_tbl[[pfs_col]] >= pfs_cut] <- "low"
        risk_group <- as.factor(risk_group)
        train_index <- create_data_partition(risk_group, p = train_prop)
    } else{
        n_train <- round(nrow(pheno_tbl) * train_prop)
        train_index <- sample(1:nrow(pheno_tbl), n_train, replace = FALSE)
    }
    # Split data by indices
    split <- list(
        "train" = list(),
        "test" = list()
    )
    split[["train"]][["expr"]] <- expr_tbl[, c(1, train_index+1)]
    split[["train"]][["pheno"]] <- pheno_tbl[train_index, ]
    split[["test"]][["expr"]] <- expr_tbl[, c(-train_index-1)]
    split[["test"]][["pheno"]] <- pheno_tbl[-train_index, ]
    if(!quiet){
        message("Splitting data into train and test data sets")
        message("Train data set has ", nrow(split[["train"]][["pheno"]]), " samples")
        message("Test data set has ", nrow(split[["test"]][["pheno"]]), " samples")
    }
    
    # Generate new DataSpec and save
    for(partition in c("train", "test")){
        ds_partition <- data_spec
        ds_partition$name <- stringr::str_c(data_spec$name, " ", partition)
        ds_partition$directory <- file.path(data_spec$directory, partition)
        split[[partition]][["data_spec"]] <- ds_partition
        if(!dir.exists(ds_partition$directory) && !quiet){
            message("Creating directory ", ds_partition$directory)
            dir.create(ds_partition$directory)
        }
        if(!quiet){
            message("Writing ", partition, " data to ", ds_partition$directory)
        }
        if(!quiet){
            message("... as ", ds_partition$expr_fname)
        }
        readr::write_csv(
            split[[partition]][["expr"]], 
            file.path(ds_partition$directory, ds_partition$expr_fname)
        )
        if(!quiet){
            message("... as ", ds_partition$pheno_fname)
        }
        readr::write_csv(
            split[[partition]][["pheno"]], 
            file.path(ds_partition$directory, ds_partition$pheno_fname)
        )
        if(!quiet){
            message("... as data_spec.rds")
        }
        saveRDS(ds_partition, file.path(ds_partition$directory, "data_spec.rds"))
    }

    invisible(split)
}