##***********************************************************************
##
##     Methods for proteasy classes
##
##***********************************************************************

#' @name proteases
#' @title Access resulting object from `findProteases` function.
#' @description `proteases` returns a `data.table` with the
#' corresponding details derived from MEROPS.
#' @aliases proteases-Cleavages proteases,Cleavages-method
#' @param x A `data.table` object.
#' @examples
#' protein <- c("P02671", "P02671", "P68871", "P01011")
#' peptide <- c("FEEVSGNVSPGTR", "FVSETESR", "LLVVYPW", "ITLLSAL")
#' res <- findProtease(protein = protein, peptide = peptide, organism = "Homo sapiens")
#' proteases(res)
#' @return A `data.table` object.
#' @rdname proteases
#' @exportMethod proteases
setMethod( f = "proteases", signature = "Cleavages", definition = function( x ) return( x@protease ))

#' @name substrates
#' @title Access resulting object from `findProteases` function.
#' @description `substrates` returns a `data.table` with the
#' corresponding details derived from MEROPS.
#' @aliases substrates-Cleavages substrates,Cleavages-method
#' @param x A `data.table` object.
#' @examples
#' protein <- c("P02671", "P02671", "P68871", "P01011")
#' peptide <- c("FEEVSGNVSPGTR", "FVSETESR", "LLVVYPW", "ITLLSAL")
#' res <- findProtease(protein = protein, peptide = peptide, organism = "Homo sapiens")
#' substrates(res)
#' @return A `data.table` object.
#' @rdname substrates
#' @exportMethod substrates
setMethod( f = "substrates", signature = "Cleavages", definition = function( x ) return( x@substrate ))

#' @name cleavages
#' @title Access resulting object from `findProteases` function.
#' @description `cleavages` returns a `data.table` with the
#' corresponding details derived from MEROPS.
#' @aliases cleavages-Cleavages cleavages,Cleavages-method
#' @param x A `data.table` object.
#' @examples
#' protein <- c("P02671", "P02671", "P68871", "P01011")
#' peptide <- c("FEEVSGNVSPGTR", "FVSETESR", "LLVVYPW", "ITLLSAL")
#' res <- findProtease(protein = protein, peptide = peptide, organism = "Homo sapiens")
#' cleavages(res)
#' @return A `data.table` object.
#' @rdname cleavages
#' @exportMethod cleavages
setMethod( f = "cleavages", signature = "Cleavages", definition = function( x ) return( x@cleavage ))
