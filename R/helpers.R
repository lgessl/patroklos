# Internal helper functions

# Convert a tibble with factor-like columns into binary indicator matrix. 
# level_list provides the possible levels for every column; the names of 
# level_list coincide with the column names of tbl.
dichotomize_tibble <- function(tbl, level_list = NULL){

    if (is.null(level_list))
        level_list <- lapply(tbl, function(c) levels(as.factor(c)))
    stopifnot(identical(colnames(tbl), names(level_list)))
    if (!all(sapply(level_list, length) > 1))
        stop("Categorical variables included from pheno data must have at ", 
            "two levels")
    # infer matrix shape and convert to factors
    n_col <- sum(sapply(level_list, function(l) length(l) - 1))
    bin_matrix <- matrix(0, nrow = nrow(tbl), ncol = n_col)

    col_idx <- 1
    bin_cnames <- character(n_col)
    for(cname in colnames(tbl)){
        lvls <- level_list[[cname]]
        for(lvl in lvls[-1]){ # first level is base level
            bin_matrix[, col_idx] <- as.numeric(tbl[[cname]] == lvl)
            # new telling name for dummy variable
            bin_cnames[col_idx] <- paste(cname, lvl, sep = "_")
            col_idx <- col_idx + 1
        }
    }
    colnames(bin_matrix) <- bin_cnames
    bin_matrix
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
    rm_na = c(FALSE, FALSE)
){
    li_var_suffix_a <- attr(a, "li_var_suffix")
    li_var_suffix_b <- attr(b, "li_var_suffix")
    if(is.matrix(a) && is.matrix(b)){
        complete_rows_a <- !logical(nrow(a))
        complete_rows_b <- !logical(nrow(b))
        if(is.null(rownames(a))) stop("intersect_by_names(): a has no row names")
        if(is.null(rownames(b))) stop("intersect_by_names(): b has no row names")
        if(rm_na[1])
            complete_rows_a <- apply(a, 1, function(row) all(!is.na(row)))
        if(rm_na[2])
            complete_rows_b <- apply(b, 1, function(row) all(!is.na(row)))
        a <- a[complete_rows_a, , drop = FALSE]
        b <- b[complete_rows_b, , drop = FALSE]
        intersect_names <- sort(intersect(rownames(a), rownames(b)))
        if(length(intersect_names) == 0){
            stop("intersect_by_names(): no common row names")
        }
        a <- a[intersect_names, , drop = FALSE]
        b <- b[intersect_names, , drop = FALSE]
        attr(a, "li_var_suffix") <- li_var_suffix_a
        attr(b, "li_var_suffix") <- li_var_suffix_b
        return(list(a, b))
    } else if(is.vector(a) && is.vector(b)){
        if(is.null(names(a))) stop("intersect_by_names(): a has no names")
        if(is.null(names(b))) stop("intersect_by_names(): b has no names")
        if(rm_na[1])
            a <- a[!is.na(a)]
        if(rm_na[2])
            b <- b[!is.na(b)]
        intersect_names <- intersect(names(a), names(b))
        if(length(intersect_names) == 0){
            stop("intersect_by_names(): no common names")
        }
        a <- a[intersect_names]
        b <- b[intersect_names]
        return(list(a, b))
    } else {
        stop("a and b must be both matrices or both vectors. However, a and b have the classes (", 
            paste(class(a), collapse = ", "), ") and (", paste(class(b), collpase = ", "), 
            "), respectively.")
    }
}

get_early_bool <- function(x, for_li = TRUE){
    li_var_suffix <- attr(x, "li_var_suffix")
    stopifnot(!is.null(colnames(x)))
    stopifnot(!is.null(li_var_suffix))
    early_bool <- vapply(
        colnames(x),
        function(s) 
            stringr::str_sub(s, -nchar(li_var_suffix)) != 
            li_var_suffix,
        logical(1)
    ) 
    if ((all(early_bool) || all(!early_bool)) && for_li) 
        stop(ifelse(all(early_bool), "All", "No"), " features are for the early", 
        " model.\n", "You set `li_var_suffix`' to ", li_var_suffix, ".\n",
        "These are your features: ", paste(colnames(x), collapse = ", "), ".\n",
        "Did you forget to set `include_from_discrete_pheno` or ", 
        "`include_from_continuous_pheno` when constructing the Model object?\n")
    early_bool
}

get_early_x <- function(x){
    early_bool <- get_early_bool(x, for_li = TRUE)
    x_early <- x[, early_bool, drop = FALSE]
    attr(x_early, "li_var_suffix") <- attr(x, "li_var_suffix")
    x_early
}

get_late_x <- function(early_predicted, x){
    early_bool <- get_early_bool(x, for_li = TRUE)
    px <- intersect_by_names(early_predicted, x[, !early_bool, drop = FALSE], 
        rm_na = c(FALSE, FALSE))
    x_late <- cbind(px[[1]], px[[2]])
    attr(x_late, "li_var_suffix") <- attr(x, "li_var_suffix")
    x_late
}