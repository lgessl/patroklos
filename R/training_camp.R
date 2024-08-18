#' @title Wrap `Model$fit()` to fit, validate and store multiple models
#' @description Given a data set and a list of `Model`s, call `Model$fit()` for 
#' all of them. If an error occurs while fitting a model, you can skip to the next 
#' model.
#' @param model_list list of `Model` objects. The models to fit.
#' @param data Data object. Its `cohort` attribute must be set. Fit the models to the `cohort` 
#' of `data`.
#' @param skip_on_error logical. Whether to skip to the next model if an error
#' occurs while fitting a model.
#' @param update_model_shell logical. If `TRUE` and, for a `Model` in `model_list`, we find a 
#' stored version with `fit_obj` not being `NULL`, we set the `fit_obj` attribute of the `Model`
#' to the found `fit_obj` and save it. This way, we can update the model shell, which we want to
#' do if changes were made to the `Model` class.
#' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
#' @export 
training_camp <- function(
    model_list,
    data,
    skip_on_error = TRUE,
    update_model_shell = FALSE,
    quiet = FALSE
){
    if(!inherits(data, "Data")){
        stop("data must be a Data object")
    }
    stopifnot(is.logical(skip_on_error))
    stopifnot(is.logical(update_model_shell))
    stopifnot(is.logical(quiet))

    if(is.null(data$expr_mat) || is.null(data$pheno_tbl)) data$read()
    if(!quiet) message("\nTRAINING CAMP ON ", data$name, ": opens at ", 
        round.POSIXt(Sys.time(), units = "secs"))
    if(is.null(data$cohort)){
        data$cohort <- "train"
        if(!quiet) message("Setting data$cohort to 'train' since it was NULL")
    }

    for(i in seq_along(model_list)){
        model <- model_list[[i]]
        if(!quiet) message("* ", model$name)
        if(!inherits(model, "Model"))
            stop("model_list must be a list of Model objects")
        if (skip_on_error) {
            tryCatch(
                {
                    model$fit(data = data, update_model_shell = update_model_shell, quiet = quiet, 
                        msg_prefix = "** ")
                }, 
                error = function(cnd){
                    warning("*** Error while fitting ", model$name, 
                        "!\n*** Error message: ", conditionMessage(cnd), 
                        "\n*** Call: ", conditionCall(cnd), 
                        "\n*** Skipping to next model.")
                }
            )
        } else {
            model$fit(data, update_model_shell = update_model_shell, quiet = quiet, 
                msg_prefix = "** ")
        }
    }
    if(!quiet) message("Training camp closes at ", round.POSIXt(Sys.time(), units = "secs"))
}