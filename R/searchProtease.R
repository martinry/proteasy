
#' @name searchProtease
#' @title Show Cleaving Data for a Peptidase or Inhibitor by Uniprot accession
#' @usage searchProtease(protein, summarize = FALSE)
#' @description Given a vector of proteins, finds which substrates they cleave.
#' @param protein a vector of corresponding UniProt Accession IDs.
#' @param summarize if false (default), provides a detailed table of all associated cleaving events,
#' otherwise outputs a summarized table and only includes reviewed (Uniprot) entries.
#'
#' @include Classes.R Generics.R Methods.R helper-functions.R
#'
#' @return data.table, character
#'
#' @examples
#' protein <- c("P98073", "P00734")
#' searchProtease(protein = protein)
#'
#' @importFrom data.table data.table
#'
#' @export

searchProtease <- function(protein, summarize = FALSE) {

    # Define local variables as NULL (due to non-standard evaluation in data.table)
    `Protease (Uniprot)` <- `Protease (MEROPS)` <- `Substrate organism` <- NULL
    `Substrate (Uniprot)` <- `Protease status` <- .N <- NULL

    # Internal data: MEROPS Substrate_search.sql and Uniprot ID to MEROPS identifier mapping
    mer <- get0("mer", envir = asNamespace("proteasy"))
    merops_map <- get0("merops_map", envir = asNamespace("proteasy"))

    # Show all data for protease

    protein_merops <- merops_map[`Protease (Uniprot)` %in% protein, c(3,4)]

    r <- mer[`Protease (MEROPS)` %in% protein_merops$`Protease (MEROPS)` & `Substrate (Uniprot)` != "\\N" &
                 `Substrate organism` %in% protein_merops$`Protease organism`]

    r <- mapMEROPSIDs(r)

    if(summarize == TRUE) r <- r[`Protease status` == "reviewed"]$`Substrate (Uniprot)`

    return(unique(r))

}
