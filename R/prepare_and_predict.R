#' @title Prepare data and predict with a loaded model across splits
#' @description Given an expression matrix and a pheno tibble, prepare the data for a
#' certain model and predict with this model. Do this for all splits into training and
#' test cohort. As with `[prepare_and_fit()]`, we do not support multiple time cutoffs
#' at this step; this is `[assessment_center()]`'s job.
#' @param expr_mat numeric matrix. The expression matrix with genes in rows and samples
#' in columns.
#' @param pheno_tbl tibble. The pheno data with samples in rows and variables in columns.
#' @param data_spec DataSpec S3 object. Specifications on the data. See the the
#' constructor `DataSpec()` for details.
#' @param model_spec ModelSpec S3 object. Specifications on the model. See the the
#' constructor `ModelSpec()` for details.
#' @param lambda string or numeric. The lambda regularization parameter of the model
#' to predict with. Technically, we will pass it to the `s` parameter of `predict()`
#' method of the object returned by the `fitter` attribute of the `ModelSpec` object.
#' See, e.g., [glmnet::predict.cv.glmnet()] or [zeroSum::predict.zeroSum()].
#' @return A list with two numeric vectors: 
#' * `"predictions"`: predicted scores,
#' *  "actual": actual, *according to `model_spec$pivot_time_cutoff` discretized* response
#' @importFrom stats predict
#' @export
prepare_and_predict <- function(
    expr_mat,
    pheno_tbl,
    data_spec,
    model_spec,
    lambda,
    benchmark = NULL
){
    if(length(model_spec$time_cutoffs) > 1L){
        stop("Multiple time cutoffs are not supported")
    }

    # Retrieve model
    fit_path <- file.path(model_spec$directory, model_spec$fit_fname)
    if(!file.exists(fit_path)){
        stop("Model object does not exist at ", fit_path)
    }
    fits <- readRDS(file.path(model_spec$directory, model_spec$fit_fname))

    model_spec$response_type <- "binary" # always evaluate for descretized response
    predicted_list <- list()
    actual_list <- list()
    benchmark_list <- list()

    for(i in model_spec$split_index){
        split_name <- paste0(data_spec$split_col_prefix, i)
        split_ms <- model_spec
        split_ms$split_index <- i
        x_y <- prepare(
            expr_mat = expr_mat,
            pheno_tbl = pheno_tbl,
            data_spec = data_spec,
            model_spec = split_ms
        )
        actual_list[[split_name]] <- x_y[["y"]][, 1]
        fit <- fits[[split_name]]
        if(is.null(fit))
            stop("No fit found for split ", split)
        predicted <- predict(fit, newx = x_y[["x"]], s = lambda)

        # Check what predict method did
        if(!is.numeric(predicted)){
            stop("predict method for class ", class(fit_obj), " does not return a ", 
            "numeric matrix or vector")
        }
        if(is.matrix(predicted)){
            if(ncol(predicted) > 1L){
                stop("predict method for class ", class(fit_obj), " returns a matrix ",
                "with more than one column")
            }
            predicted <- predicted[, 1]
        }
        if(is.null(names(predicted))){
            names(predicted) <- rownames(x_y[["x"]])
        }
        predicted_list[[split_name]] <- predicted
        benchmark_list[[split_name]] <- pheno_tbl[[benchmark]][names(predicted)]
    }

    res <- list(
        "predicted" = predicted_list, 
        "actual" = actual_list,
        "benchmark" = benchmark_list
    )
    return(res)
}