#' @export 
training_camp <- function(
    data_spec,
    model_spec_list
){
    # Read in data once and for all
    if(!inherits(data_spec, "DataSpec")){
        stop("data_spec must be a DataSpec object")
    }
    data <- read(data_spec)
    expr_mat <- data[["expr_mat"]]
    pheno_tbl <- data[["pheno_tbl"]]

    pheno_tbl <- ensure_splits(
        pheno_tbl = pheno_tbl,
        data_spec = data_spec,
        model_spec_list = model_spec_list
    )

    for(i in seq_along(model_spec_list)){

        model_spec <- model_spec_list[[i]]
        if(!inherits(model_spec, "ModelSpec")){
            stop("model_spec_list must be a list of ModelSpec objects")
        }
        for(j in seq_along(model_spec$cutoff_times)){
            # Adapt ModelSpec to very cutoff_time
            ms_cutoff <- model_spec # Rely on copy-on-modify
            cutoff_time <- model_spec$cutoff_times[j]
            ms_cutoff$name <- paste0(model_spec$name, "@", cutoff_time)
            ms_cutoff$cutoff_times <- cutoff_time
            ms_cutoff$directory <- file.path(
                model_spec$directory, 
                stringr::str_replace(as.character(cutoff_time), "\\.", "-")
            )
            prepare_and_fit(
                expr_mat = expr_mat,
                pheno_tbl = pheno_tbl,
                data_spec = data_spec,
                model_spec = ms_cutoff
            )
        }
    }
}