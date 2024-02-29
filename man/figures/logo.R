p <- ggplot2::ggplot() + ggplot2::theme_void()

hexSticker::sticker(p, package = "patroklos", h_fill = "black", h_color = "white", 
    p_size = 23, p_x = 1, p_y = 1.1, filename = "man/figs/logo_raw.png", 
    spotlight = FALSE, dpi = 300)
