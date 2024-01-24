#' @title Set up a training camp to fit models
#' @description Given one data set and a list of models, fit all models for all 
#' splits and time cutoffs. Store the models.
#' @param data_spec DataSpec object. The data set to fit on.
#' @param model_spec_list list of ModelSpec objects. The models to fit.
#' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
#' @return NULL
#' @export 
training_camp <- function(
    data_spec,
    model_spec_list,
    quiet = FALSE
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

    if(!quiet) message("\nTRAINING CAMP ON ", data_spec$name, ": opens at ", 
        round(Sys.time(), units = "secs"))
    for(i in seq_along(model_spec_list)){
        model_spec <- model_spec_list[[i]]
        if(!quiet) message("# ", model_spec$name)
        if(!inherits(model_spec, "ModelSpec")){
            stop("model_spec_list must be a list of ModelSpec objects")
        }
        for(time_cutoff in model_spec$time_cutoffs){
            if(!quiet) message("## At time cutoff ", time_cutoff, 
                " (", round(Sys.time(), units = "secs"), ")")
            ms_cutoff <- at_time_cutoff(model_spec, time_cutoff)
            prepare_and_fit(
                expr_mat = expr_mat,
                pheno_tbl = pheno_tbl,
                data_spec = data_spec,
                model_spec = ms_cutoff,
                quiet = quiet,
                msg_prefix = "### "
            )
        }
    }
    if(!quiet) message("Training camp closes at ", round(Sys.time(), units = "secs"))
}