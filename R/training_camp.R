#' @title Set up a training camp to fit models
#' @description Given one data set and a list of models, fit all models for all 
#' splits and time cutoffs. Store the models. If an error occurs while fitting 
#' a model, skip to next model.
#' @param model_list list of Model objects. The models to fit.
#' @param data Data object. The data set to fit on.
#' @param skip_on_error logical. Whether to skip to the next model if an error
#' occurs while fitting a model.
#' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
#' @export 
training_camp <- function(
    model_list,
    data,
    skip_on_error = TRUE,
    quiet = FALSE
){
    if(!inherits(data, "Data")){
        stop("data must be a Data object")
    }
    stopifnot(is.logical(quiet))

    if(is.null(data$expr_mat) || is.null(data$pheno_tbl)) data$read()
    if(!quiet) message("\nTRAINING CAMP ON ", data$name, ": opens at ", 
        round.POSIXt(Sys.time(), units = "secs"))
    if(is.null(data$cohort)){
        data$cohort <- "train"
        if(!quiet) message("Setting data$cohort to 'train' since it was NULL")
    }
    data$pheno_tbl <- ensure_splits(
        pheno_tbl = data$pheno_tbl,
        data = data,
        model_list = model_list
    )

    for(i in seq_along(model_list)){
        model <- model_list[[i]]
        if(!quiet) message("* ", model$name)
        if(!inherits(model, "Model"))
            stop("model_list must be a list of Model objects")
        if (skip_on_error) {
            tryCatch(
                {
                    model$fit(data, quiet = quiet, msg_prefix = "** ")
                }, 
                error = function(cnd){
                    warning("*** Error while fitting ", model$name, 
                        "!\n*** Error message: ", conditionMessage(cnd), 
                        "\n*** Call: ", conditionCall(cnd), 
                        "\n*** Skipping to next model.")
                }
            )
        } else {
            model$fit(data, quiet = quiet, msg_prefix = "*** ")
        }
    }
    if(!quiet) message("Training camp closes at ", round.POSIXt(Sys.time(), units = "secs"))
}