# Internal helper functions

# Convert a tibble with factor-like columns into binary indicator matrix 
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


# Are all elements of a vector unique? FALSE if NA is present.
elements_unique <- function(x){
    if(any(is.na(x))) return(FALSE)
    return(length(unique(x)) == length(x))
}


# Split the indices of the factor `outcome` into two groups by sampling a 
# fraction `p` of the indices belonging to each level of `outcome`.
create_data_partition <- function(
    outcome,
    p
){
    if(!is.factor(outcome)){
        stop("outcome must be a factor")
    }

    selected <- NULL
    for(lvl in levels(outcome)){
        at_stake <- which(outcome == lvl)
        selected <- c(
            selected, 
            sample(at_stake, round(p*length(at_stake)), replace = FALSE)
        )
    }
    return(selected)
}


# For two matrices or two vectors, intersect the (row)names and return the
# corresponding ordered subset of both matrices or vectors. Additonally, you can
# remove values (rows) with NA 
intersect_by_names <- function(
    a,
    b,
    rm_na = FALSE
){
    if(is.matrix(a) && is.matrix(b)){
        if(is.null(rownames(a))) stop("intersect_by_names(): a has no row names")
        if(is.null(rownames(b))) stop("intersect_by_names(): b has no row names")
        if(rm_na){
            complete_rows_a <- apply(a, 1, function(row) all(!is.na(row)))
            complete_rows_b <- apply(b, 1, function(row) all(!is.na(row)))
            a <- a[complete_rows_a, , drop = FALSE]
            b <- b[complete_rows_b, , drop = FALSE]
        }
        intersect_names <- intersect(rownames(a), rownames(b)) |> sort()
        if(length(intersect_names) == 0){
            stop("intersect_by_names(): no common row names")
        }
        a <- a[intersect_names, , drop = FALSE]
        b <- b[intersect_names, , drop = FALSE]
        return(list(a, b))
    } else if(is.vector(a) && is.vector(b)){
        if(is.null(names(a))) stop("intersect_by_names(): a has no names")
        if(is.null(names(b))) stop("intersect_by_names(): b has no names")
        if(rm_na){
            a <- a[!is.na(a)]
            b <- b[!is.na(b)]
        }
        intersect_names <- intersect(names(a), names(b))
        if(length(intersect_names) == 0){
            stop("intersect_by_names(): no common names")
        }
        a <- a[intersect_names]
        b <- b[intersect_names]
        return(list(a, b))
    } else {
        stop("intersect_by_names(): a and b must be both matrices or both vectors")
    }
}


# Replace unique directory in filepath with another directory
mirror_directory <- function(
    filepath,
    mirror
){
    path_split <- stringr::str_split(filepath, .Platform$file.sep)[[1]]
    replace_at <- which(path_split == mirror[1])
    if(length(replace_at) != 1)
        stop("mirror_directory(): mirror[1] must match exactly one element of filepath",
            "\n filepath: ", filepath, "\n mirror: ", mirror)
    path_split[replace_at] <- mirror[2]
    new_filepath <- do.call(file.path, as.list(path_split))
    return(new_filepath)
}