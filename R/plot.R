#' @importFrom rlang .data
plot_2d_metric <- function(
    tbl,
    ass2d,
    quiet = FALSE,
    msg_prefix = ""
){
    plt <- ggplot2::ggplot(
        data = tbl,
        mapping = ggplot2::aes(
            x = .data[[ass2d$x_metric]], 
            y = .data[[ass2d$y_metric]]
        )
    ) + 
        ass2d$theme +
        ggplot2::geom_point(alpha = ass2d$alpha, shape = 4) +
        ggplot2::labs(
            title = ass2d$title, 
            x = ass2d$x_lab, 
            y = ass2d$y_lab
        ) +
        ggplot2::scale_x_continuous(trans = ass2d$scale_x) + 
        ggplot2::scale_y_continuous(trans = ass2d$scale_y)
    if (length(unique(tbl[["model"]])) > 1) plt <- plt + ggplot2::aes(color = .data[["model"]])
    if(!is.null(ass2d$colors)){
        plt <- plt + ggplot2::scale_color_manual(values = ass2d$colors)
    }

    if(!is.null(ass2d$file)){
        if(!quiet)
            message(msg_prefix, "Saving 2D metric plot to ", ass2d$file)
        showtext::showtext_auto()
        ggplot2::ggsave(
            ass2d$file,
            plt,
            width = ass2d$width, 
            height = ass2d$height, 
            units = ass2d$units,
            dpi = ass2d$dpi
        )
        showtext::showtext_auto(FALSE)

        # Save to csv (if wanted)
        if(ass2d$fellow_csv){
            csv_file <- stringr::str_replace(ass2d$file, "\\..+", ".csv")
            if(!quiet) message(msg_prefix, "Saving 2D metric table to ", csv_file)
            readr::write_csv(tbl, csv_file)
        }
    }
    return(plt)
}