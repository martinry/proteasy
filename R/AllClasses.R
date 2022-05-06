
setClass(
	Class = "Cleavages",
	representation = list(
		organism   = "character",
		substrate  = "data.table",
		protease   = "data.table",
		cleavage   = "data.table"
	)
)


setGeneric("proteases", function(x) standardGeneric("proteases"))
setGeneric("substrates", function(x) standardGeneric("substrates"))
setGeneric("cleavages", function(x) standardGeneric("cleavages"))
setGeneric("browseProtease", function(x,y) standardGeneric("browseProtease"))

`.` <- list
