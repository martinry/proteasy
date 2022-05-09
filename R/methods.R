##***********************************************************************
##
##     Methods for proteasy classes
##
##***********************************************************************

#' @title Access resulting object from `findProteases` function.
#' @description `proteases`, `substrates` and `cleavages` returns a `data.table` with the
#' corresponding details derived from MEROPS.
#'
#' @param object A `data.table` object.
#'
#' @return A `data.table` object.
#'
#' @export
#' @md

setMethod( f = "proteases", signature = "Cleavages", definition = function( x ) return( x@protease ))

#' @title Access resulting object from `findProteases` function.
#' @description `proteases`, `substrates` and `cleavages` returns a `data.table` with the
#' corresponding details derived from MEROPS.
#'
#' @param object A `data.table` object.
#'
#' @return A `data.table` object.
#'
#' @export
#' @md

setMethod( f = "substrates", signature = "Cleavages", definition = function( x ) return( x@substrate ))

#' @title Access resulting object from `findProteases` function.
#' @description `proteases`, `substrates` and `cleavages` returns a `data.table` with the
#' corresponding details derived from MEROPS.
#'
#' @param object A `data.table` object.
#'
#' @return A `data.table` object.
#'
#' @export
#' @md

setMethod( f = "cleavages", signature = "Cleavages", definition = function( x ) return( x@cleavage ))
