import::here(stringi, 'stri_replace_all_regex')

## Functions
## chunker
## collate
## dict_zip
## filter_list_for_match
## index_last_occurrence
## items_in_a_not_b
## multiple_replacement
## replace_specific_items


#' Chunker
#' 
#' @description Group array elements into chunks for iteration
#' @examples
#' chunker(c(1, 2, 3, 4, 5, 6), 4)
#' list(c(1, 2, 3, 4), c(5, 6))
#' 
chunker <- function(array_obj, num_elem=2) {

    num_groups = ceiling(length(array_obj)/num_elem)
    chunked <- vector("list", length = num_groups)
    for (group_idx in 1:num_groups) {
        start <- 1+(group_idx-1)*num_elem
        stop <- ifelse(group_idx*num_elem > length(array_obj), length(array_obj), group_idx*num_elem)
        chunked[[group_idx]] <- array_obj[start:stop]
    }
    return(chunked)
}


#' Collate
#' 
#' @description Collate array elements for iteration
#' @examples
#' collate(c(1, 2, 3, 4, 5, 6), 4)
#' list(c(1, 5), c(2, 6), 3, 4)
#' 
collate <- function(array_obj, num_groups=2) {

    if (num_groups >= length(array_obj)) {
        return(as.list(array_obj))
    }

    collated <- vector("list", length = num_groups)
    for (group_idx in 1:num_groups) {
        collated[[group_idx]] <- array_obj[seq(group_idx, length(array_obj), num_groups)]
    }
    return(collated)
}


#' Dictionary
#' 
#' @description Simple dictionary implementation using R environment
#'
#' @export
dict_zip <- function(keys, values) {
    if (!is.list(values)) {
        values <- as.list(values)
    }
    named_list <- setNames(values, keys)
    env <- list2env(named_list)
    return(env)
}


#' Find matches based on substring
#' 
#' @description Return elements of a list matching a particular substring
#' 
#' @param items list or vector
#' @param patterns a string or collection of strings
#' @return Returns a list or vector with any matching items
#' 
#' @examples
#' filter_list_for_match(c("a_suffix", "b_suffix", "c"), "suffix")
#' 
filter_list_for_match <- function(items, pattern) {
    # filter
    for (i in 1:length(pattern)){
        items <- lapply(items, grep, pattern=pattern[[i]], value=TRUE)
    }
    return (unlist(items[!sapply(items, identical, character(0))]))  # remove character(0)
}


#' Index Last Occurrence
#' 
#' @description
#' Returns a named list of the last occurrence of an item in a list
#' 
#' @examples
#' index_last_occurrence(c('a', 'a', 'a', 'b', 'b', 'c'), c('a', 'b', 'c'))
#' 
index_last_occurrence <- function(items) {
    indexes <- sort(tapply(seq_along(items), items, max))
    return(indexes)
}


#' Return items unique to vector a
#' 
#' @description
#' Return all the items found in a and not b. Useful for filtering dataframe columns.
#' 
#' @param a list or vector
#' @param b list or vector
#' @return Returns the filtered list or vector (same type as a)
#' 
#' @references
#' \href{https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list}{Stack Overflow}
#' 
#' @examples
#' items_in_a_not_b(c('a', 'b', 'c', '1', '2'), c('1', '2'))
#' 
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}


#' Replaces each item in a list or vector with all replacements
#' 
#' @description
#' Convenience function to perform multiple replacements on a list or dataframe column.
#' Unlike [replace_specific_items()], `multiple_replacement()` can recognize patterns.
#' 
#' @param items list or vector
#' @param replacements a named list of replacements. uses names to match and values to replace.
#' @return Returns a vector with replaced items
#' 
#' @examples
#' replacements <- c('prefix_' = '', '_suffix' = '')
#' items <- c('prefix_a_suffix', 'prefix_b_suffix')
#' multiple_replacement(items, replacements)
#' 
#' @seealso [replace_specific_items()], [stringi::stri_replace_all_regex()]
#' 
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


#' Replaces each item in a list or vector with all replacements
#' 
#' @description Use this to rename columns
#' 
#' @param items list or vector
#' @param replacements a named list of replacements. uses names to match and values to replace.
#' @return Returns a vector with replaced items
#' 
#' @examples
#' replace_specific_items(c('a', 'b', 'c'), c('a'="A", 'c'="C"))
#' 
#' @seealso [multiple_replacement()]
#' 
#' @export
replace_specific_items <- function(items, replacements) {
    replace_ids <- which(items %in% intersect(names(replacements), items))
    for (idx in replace_ids) {
        items[idx] <- replacements[items[idx]]
    }
    return(items)
}
