# convert a tibble with factor-like columns into binary indicator matrix 
# using dummy variables
tibble_to_binary <- function(
    tbl
){
    # type checking
    if(!is.data.frame(tbl)){
        stop("tibble_to_binary: tbl must be a tibble or data frame")
    }

    # infer matrix shape and convert to factors
    ncol <- 0
    for(cname in colnames(tbl)){
        fac <- as.factor(tbl[[cname]])
        tbl[[cname]] <- fac
        ncol <- ncol + length(levels(fac)) - 1
    }
    bin_matrix <- matrix(0, nrow = nrow(tbl), ncol = ncol)

    col_idx <- 1
    bin_cnames <- character(ncol)
    for(cname in colnames(tbl)){
        lvls <- levels(tbl[[cname]])
        for(lvl in lvls[2:length(lvls)]){ # first level is base level
            bin_matrix[, col_idx] <- as.numeric(tbl[[cname]] == lvl)
            # new telling name for dummy variable
            bin_cnames[col_idx] <- paste(cname, lvl, sep = "_")
            col_idx <- col_idx + 1
        }
    }
    colnames(bin_matrix) <- bin_cnames
    return(bin_matrix)
}


elements_unique <- function(x){
    return(length(unique(x)) == length(x))
}