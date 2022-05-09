##***********************************************************************
##
##     proteasy classes
##
##***********************************************************************

setClass(
    Class = "Cleavages",
    slots = c(
        organism   = "character",
        substrate  = "data.table",
        protease   = "data.table",
        cleavage   = "data.table"
    )
)
