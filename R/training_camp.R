#' @title Set up a training camp to fit models
#' @description Given one data set and a list of models, fit all models for all 
#' splits and time cutoffs. Store the models.
#' @param model_list list of ModelSpec objects. The models to fit.
#' @param data DataSpec object. The data set to fit on.
#' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
#' @return NULL
#' @export 
training_camp <- function(
    model_list,
    data,
    quiet = FALSE
){
    # Read in data once and for all
    if(!inherits(data, "DataSpec")){
        stop("data must be a DataSpec object")
    }
    data <- read(data)
    expr_mat <- data[["expr_mat"]]
    pheno_tbl <- data[["pheno_tbl"]]

    pheno_tbl <- ensure_splits(
        pheno_tbl = pheno_tbl,
        data = data,
        model_list = model_list
    )

    if(!quiet) message("\nTRAINING CAMP ON ", data$name, ": opens at ", 
        round.POSIXt(Sys.time(), units = "secs"))
    for(i in seq_along(model_list)){
        model <- model_list[[i]]
        if(!quiet) message("# ", model$name)
        if(!inherits(model, "ModelSpec")){
            stop("model_list must be a list of ModelSpec objects")
        }
        for(time_cutoff in model$time_cutoffs){
            if(!quiet) message("## At time cutoff ", time_cutoff, 
                " (", round.POSIXt(Sys.time(), units = "secs"), ")")
            model_cutoff <- at_time_cutoff(model, time_cutoff)
            prepare_and_fit(
                expr_mat = expr_mat,
                pheno_tbl = pheno_tbl,
                data = data,
                model = model_cutoff,
                quiet = quiet,
                msg_prefix = "### "
            )
        }
    }
    if(!quiet) message("Training camp closes at ", round.POSIXt(Sys.time(), units = "secs"))
}