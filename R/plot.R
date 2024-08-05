#' @importFrom rlang .data
plot_2d_metric <- function(
    tbl,
    ass2d,
    file,
    fellow_csv,
    quiet = FALSE,
    msg_prefix = ""
){
    plt <- ggplot2::ggplot(
        data = tbl,
        mapping = ggplot2::aes(
            x = .data[[ass2d$x_metric]], 
            y = .data[[ass2d$y_metric]],
            color = .data[["model"]]
        )
    ) + 
        ass2d$theme +
        ggplot2::labs(
            x = ass2d$x_lab, 
            y = ass2d$y_lab
        ) +
        ggplot2::scale_x_continuous(trans = ass2d$scale_x) + 
        ggplot2::scale_y_continuous(trans = ass2d$scale_y) +
        ggplot2::geom_line(alpha = ass2d$alpha)

    if(!is.null(file)){
        if(!quiet)
            message(msg_prefix, "Saving 2D metric plot to ", file)
        showtext::showtext_auto()
        ggplot2::ggsave(
            file,
            plt,
            width = ass2d$width, 
            height = ass2d$height, 
            units = ass2d$units,
            dpi = ass2d$dpi
        )
        showtext::showtext_auto(FALSE)

        # Save to csv (if wanted)
        if(fellow_csv){
            csv_file <- stringr::str_replace(file, "\\..+", ".csv")
            if(!quiet) message(msg_prefix, "Saving 2D metric table to ", csv_file)
            readr::write_csv(tbl, csv_file)
        }
    }
    return(plt)
}