qc_prepare <- function(
    x,
    y
){
    # check class and type of x and y
    x_y_list <- list("x" = x, "y" = y)
    for(x_y in names(x_y_list)){
        if(!("matrix" %in% class(x_y_list[[x_y]]))){
            stop(x_y, "must be a matrix.")
        }
        if(typeof(x_y_list[[x_y]]) !="double"){
            stop(x_y, "must be of type double.")
        }
    }

    check_consistent_patient_ids(
        stage = "after_generate_xy",
        expr = x,
        pheno = y
    )
    check_available(
        x = x,
        y = y
    )
}