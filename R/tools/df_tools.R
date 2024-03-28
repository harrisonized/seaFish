import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'replace_specific_items', .character_only=TRUE)

## Functions
## fillna
## matrix_to_dataframe
## smelt
## set_index
## rename_columns
## reset_index
## rev_df
## dgCMatrix_to_dataframe


#' Fill specific column with NA
#' 
#' @description Mirrors Pandas' \href{https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.fillna.html}{fillna}
#' 
#' @param df a dataframe
#' @param cols a list of columns to replace
#' @param val the value to fill with
#' @param inplace TRUE allows you to avoid re-assigning the variable
#' @return Returns a dataframe.
#' 
#' @examples
#' mtcars['new_col'] <- NA
#' head(fillna(mtcars, c('new_col'), 1))
#' 
#' @export
fillna <- function(df, cols, val=0, inplace=FALSE) {
    df_name <- deparse(substitute(df))
    for (col in cols) {
        df[is.na(df[, col]), col] <- val
    }
    if (inplace) {
        assign(df_name, df, envir=.GlobalEnv)
    } else {
        return(df)
    }
}


#' Matrix to Dataframe
#'
matrix_to_dataframe <- function(
    matrix
) {
    return (as.data.frame(as.matrix(expr_mtx)))
}


#' Special Melt
#' 
#' @description
#' Convert a dataframe from a table to long format
#' This is an alternative to melt that doesn't throw errors
#' 
#' @examples
#' smelt(mtcars[, c('cyl', 'mpg')])
#' 
#' @references
#' See: \href{https://stackoverflow.com/questions/28355377/how-to-add-index-of-a-list-item-after-melt-in-r}{Stack Overflow link}
#' 
#' @export
smelt <- function(
   df,
   rowname='row',
   colname='col',
   valname='val'
) {
   melted <- transform(stack(setNames(df, colnames(df))), id=rownames(df))
   colnames(melted) <- c(valname, colname, rowname)
   return(rev(melted))
}


#' Set Index
#' 
#' @description
#' Uses the values in a column to set an index. Needs to be unique.
#' Mirrors Pandas' \href{https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.set_index.html}{set_index}.
#' 
#' @export
set_index <- function(df, colname, drop=TRUE) {

    rownames(df) <- df[[colname]]
    if (drop) {
        df <- df[, !names(df)==colname]
    }
    return(df)
}


#' Rename dataframe columns
#' 
#' @description This is a thin wrapper around replace_specific_items that acts on dataframe columns
#' 
#' @param df a dataframe
#' @param columns a named list of replacements. uses names to match and values to replace
#' @param inplace TRUE allows you to avoid re-assigning the variable
#' @return Returns a dataframe.
#' 
#' @examples
#' head(rename_columns(mtcars, c('mpg'="MPG", 'disp'="DISP")))
#' 
#' @export
rename_columns <- function(df, columns, inplace=FALSE) {
    df_name <- deparse(substitute(df))
    colnames(df) <- replace_specific_items(colnames(df), columns)
    if (inplace) {
        assign(df_name, df, envir=.GlobalEnv)
    } else {
        return(df)
    }
}


#' Reset index
#' 
#' @description
#' Moves the values in the index to a column. Resets the index to the default integer index.
#' Mirrors Pandas' \href{https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.reset_index.html}{reset_index}.
#' 
#' @param df a dataframe
#' @param index_name select a new name for the index column
#' @param drop if TRUE, does not copy the index values to the new column
#' @return Returns a dataframe with index renamed and new integer index 
#' 
#' @examples
#' reset_index(mtcars, index_name='model')
#' 
#' @references
#' \href{https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column}{StackOverflow post}
#' 
#' @export
reset_index <- function(df, index_name='index', drop=FALSE) {
    if (!drop) {
        df <- cbind(index = rownames(df), df)
        colnames(df)[colnames(df) == "index"] = index_name
    }
    rownames(df) <- 1:nrow(df)
    return (df)
}


#' Reverse the order of a dataframe
#' 
#' @param df a dataframe
#' @param how 'row' or 'col'
#' @return Returns the reversed dataframe.
#' 
#' @examples
#' rev_df(mtcars, how='row')
#' rev_df(mtcars, how='col')
#' 
#' @export
rev_df <- function(df, how='row') {
    if (how == 'row') {
        return(df[dim(df)[1]:1,])
    } else if (how == 'col') {
        return(rev(df))
    } else {
        stop("Choose how='row' or how='col'")
    }
}
