import::here(digest, 'sha1')

## Functions
## title_to_snake_case
## random_hash


#' Random Hash
#' 
#' @description Generates a Github style hash
#' 
#' @examples
#' random_hash()
#' 
#' @export
random_hash <- function(digits=6) {
    return(substr(sha1(runif(1, 1, 2^31-1), digits = 14), 1, digits))
}


#' Standardize Space-separated Titles
#' 
#' @description Converts "Column Title" to column_title
#' 
#' @examples
#' title_to_snake_case('Column Title')
#' 
#' @export
title_to_snake_case <- function(text) {
    return(tolower(
        paste(
            unlist(strsplit(text, '[ ]')), collapse='_')
        )
    )
}
