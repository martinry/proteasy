

setMethod( f = "proteases", signature = "Cleavages", definition = function( x ) return( x@protease ))
setMethod( f = "substrates", signature = "Cleavages", definition = function( x ) return( x@substrate ))
setMethod( f = "cleavages", signature = "Cleavages", definition = function( x ) return( x@cleavage ))
