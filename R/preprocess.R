# Helpful functions for preprocessing

#' @title Add discretized IPI features to a pheno_tbl
#' @description Discretize the columns in a pheno tibble that belong to the features of
#' the International Prognostic Index (IPI) into binary indicator columns by comparing them 
#' to cutoff values and append them to the pheno tibble. 
#' @param pheno_tbl tibble. The pheno data.
#' @param data_spec DataSpec. The data specification. See `DataSpec` for details.
#' @param cutoffs named numeric vector. The cutoff values for the IPI features. The names
#' of it must coincide with the names of `data_spec$ipi_feat_cols`. The default (see Usage)
#' takes the cutoffs from the original IPI paper [
#' The International Non-Hodgkin's Lymphoma Prognostic Factors Project (1993):
#' A Predictive Model for Aggressive Non-Hodgkin's Lymphoma
#' ](https://www.nejm.org/doi/full/10.1056/NEJM199309303291402).
#' @return The pheno tibble with the discretized IPI features appended.
#' @export
discretize_ipi_features <- function(
    pheno_tbl,
    data_spec,
    cutoffs = c(
        "age" = 60.,
        "stage" = 2.,
        "ldh_ratio" = 1.,
        "performance_status" = 1.,
        "n_extranodal_sites" = 1.
    )
){
    cutoff_names <- names(cutoffs)
    ipi_feat_names <- names(data_spec$ipi_feat_cols)
    # Checks
    if(!setequal(cutoff_names, ipi_feat_names)){
        stop("discretize_ipi_features():  cutoffs must have the same names as ",
            "`data_spec$ipi_feat_cols`, namely ", 
            paste0(ipi_feat_names, collapse = ", "))
    }
    # Extend pheno_tbl
    for(feat in ipi_feat_names){
        feat_colname <- data_spec$ipi_feat_cols[feat]
        feat_col <- pheno_tbl[[feat_colname]]
        if(is.null(feat_col)){
            stop("discretize_ipi_features():  pheno_tbl has no column named ",
                feat_colname)
        }
        pheno_tbl[[paste0(feat_colname, ">", cutoffs[feat])]] <- as.numeric(
            feat_col > cutoffs[feat]
        )
    }
    return(pheno_tbl)
}