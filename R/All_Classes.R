

setClass("proteasy_res", representation(
	organism   = "character",
	substrate  = "data.table",
	protease   = "data.table",
    output     = "data.table"
    )
)