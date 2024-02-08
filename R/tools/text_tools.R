import::here(digest, 'sha1')

## Functions
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
