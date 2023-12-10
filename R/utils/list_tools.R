import::here('stringi', 'stri_replace_all_regex', .character_only=TRUE)

## Functions
## items_in_a_not_b
## filter_list_for_match
## multiple_replacement


#' https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list
#' 
#' @export
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}


#' return elements of a list matching a particular substring
#'
#' @examples
#' filter_list_for_match(c("gene_id_pat", "gene_id_mat", "count"), "pat")
#' 
#' @export
filter_list_for_match <- function(items, pattern) {
    # filter
    for (i in 1:length(pattern)){
        items <- lapply(items, grep, pattern=pattern[[i]], value=TRUE)
    }
    return (unlist(items[!sapply(items, identical, character(0))]))  # remove character(0)
}


#' Convenience function to perform multiple replacements on a list or dataframe column
#' 
#' Example:
#' replace_dict <- c(
#'     '[A-Za-z]' = '',
#'     '[0-9]+' = ''
#' )
#' 
#' @export
multiple_replacement <- function(items, replace_dict) {

    patterns <- names(replace_dict)
    replacements <- sapply(unname(replace_dict), function(x) gsub('\\\\', '$', x))
    
    items <- sapply(items,
        function(x) stri_replace_all_regex(
            x,
            pattern = patterns,
            replacement = replacements,
            vectorize_all = FALSE)
    )

    return (items)
}
