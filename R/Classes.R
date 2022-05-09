##***********************************************************************
##
##     proteasy classes
##
##***********************************************************************

#' @slot organism a character vector describing the organism of substrates and proteases in the object.
#' @slot substrate a data.table detailing substrate data.
#' @slot protease a data.table detailing protease data.
#' @slot cleavage a data.table detailing cleavage data.
setClass(
    Class = "Cleavages",
    slots = c(
        organism   = "character",
        substrate  = "data.table",
        protease   = "data.table",
        cleavage   = "data.table"
    )
)
