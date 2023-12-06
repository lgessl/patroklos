#' @title Download and read data into a tibble
#' @description Download a csv, tsv or xlsx file from a url (if not already 
#' downloaded), read it in and return it as a tibble.
#' @param data_dir string. The directory to store the downloaded file in.
#' @param url string. The url to download the file from.
#' @param to_file string. Dowloaded file will be stored as `to_file` inside
#' `data_dir`.
#' @param ... Additional arguments to pass to the reader function. For example,
#' `readr::read_csv()` takes `col_types` as an argument.
#' @return A tibble.
#' @export
download_and_read <- function(
    data_dir,
    url,
    to_file,
    ...
){
    # Download
    if(!dir.exists(data_dir)){
        message("Creating directory", data_dir)
        dir.create(data_dir)
    }
    to_file <- file.path(data_dir, to_file)
    if(!file.exists(to_file)){
        message("Downloading", url)
        download.file(url, to_file)
    }
    # Read
    file_extension <- stringr::str_extract(to_file, "\\..+$")
    reader <- switch(
        stringr::str_extract(to_file, "\\..+$"),
        ".csv" = readr::read_csv,
        ".tsv" = readr::read_tsv,
        ".xlsx" = readxl::read_excel,
        stop("Only support csv, tsv, and xlsx files.")
    )
    message("Reading in ", to_file)
    res <- do.call(
        reader, 
        args = c(list(to_file), list(...))
    )
    return(res)
}