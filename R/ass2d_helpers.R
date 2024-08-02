ass2d_initialize <- function(self, private, x_metric, y_metric, 
    ci_level, x_lab, y_lab, xlim, ylim, scale_x, scale_y, 
    alpha, colors, theme, width, height, units, dpi)
{
    if(is.null(x_lab)){
        x_lab <- x_metric
    }
    if(is.null(y_lab)){
        y_lab <- y_metric
    }
    if(y_metric == "logrank"){
        if(!any(stringr::str_detect(c("prevalence", "rpp"), x_metric)))
            stop("For `y_metric` = 'logrank', `x_metric` must be 'prevalence' or 'rpp'.")
    }
    stopifnot(is.character(x_metric))
    stopifnot(is.character(y_metric))
    stopifnot(is.numeric(ci_level) && ci_level >= 0 && ci_level <= 1)
    stopifnot(is.character(x_lab))
    stopifnot(is.character(y_lab))
    stopifnot(is.numeric(xlim) || is.null(xlim) || 
        xlim[1] >= xlim[2] || ylim[1] >= ylim[2])
    stopifnot(is.numeric(ylim) || is.null(ylim))
    stopifnot(is.numeric(alpha) && alpha >= 0 && alpha <= 1)
    stopifnot(is.character(colors) || is.null(colors))
    stopifnot(is.numeric(width) && width > 0)
    stopifnot(is.numeric(height) && height > 0)
    stopifnot((inherits(theme, "theme") && inherits(theme, "gg")) || is.null(theme))
    stopifnot(is.character(units))
    stopifnot(is.numeric(dpi) && dpi > 0)

    self$x_metric <- x_metric
    self$y_metric <- y_metric
    self$ci_level <- ci_level
    self$x_lab <- x_lab
    self$y_lab <- y_lab
    self$xlim <- xlim
    self$ylim <- ylim
    self$scale_x <- scale_x
    self$scale_y <- scale_y
    self$alpha <- alpha
    self$colors <- colors
    self$theme <- theme
    self$width <- width
    self$height <- height
    self$units <- units
    self$dpi <- dpi

    invisible(self)
}

ass2d_assess <- function(self, private, data, model, return_type, file, fellow_csv, quiet, 
    msg_prefix){

    stopifnot(inherits(data, "Data"))
    stopifnot(inherits(model, "Model"))
    stopifnot(is.character(return_type))
    stopifnot(is.null(file) || is.character(file))
    stopifnot(is.logical(fellow_csv))
    stopifnot(is.logical(quiet))
    stopifnot(is.character(msg_prefix))
    return_type <- match.arg(return_type, c("ggplot", "tbl_df"))
    if (!is.null(file)) {
        directory <- dirname(file)
        if(!dir.exists(directory))
            dir.create(directory, recursive = TRUE)
    }
    if(is.null(data$expr_mat) || is.null(data$pheno_tbl)) data$read()
    tbl <- private$get_2d_metric(
        data = data, 
        model = model,
        quiet = quiet
    )
    if (return_type == "tbl_df") return(tbl)
    if (is.null(tbl)) return(NULL)
    plt <- plot_2d_metric(
        tbl = tbl,
        ass2d = self,
        file = file,
        fellow_csv = fellow_csv,
        msg_prefix = msg_prefix,
        quiet = quiet
    )
    return(plt)
}

ass2d_assess_center <- function(self, private, data, model_list, file, fellow_csv, quiet){
    
    stopifnot(inherits(data, "Data"))   
    stopifnot(is.list(model_list))  
    stopifnot(all(sapply(model_list, function(x) inherits(x, "Model"))))
    stopifnot(is.logical(quiet))
    if (is.null(data$expr_mat) || is.null(data$pheno_tbl)) data$read()
    if(!quiet) message("Assessment center on ", data$name, " open")
    tbls <- vector("list", length(model_list))

    # Assess each model
    for (i in seq_along(model_list)) {
        model <- model_list[[i]]
        if (!is.null(model$continuous_output) && !model$continuous_output)
            next
        if(!quiet) message("** ", model$name)
        tbls[[i]] <- self$assess(
            data = data,
            model = model,
            return_type = "tbl_df",
            quiet = quiet,
            msg_prefix = "**** "
        )
    }
    tbls <- tbls[!sapply(tbls, is.null)]
    tbl <- dplyr::bind_rows(tbls)

    plt <- plot_2d_metric(
        tbl = tbl,
        ass2d = self,
        file = file,
        fellow_csv = fellow_csv,
        quiet = quiet
    )
    return(plt)
}