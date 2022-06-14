
#' @name searchSubstrate
#' @title Show Cleaving Data for a Substrate by Uniprot accession
#' @usage searchSubstrate(protein, summarize = FALSE)
#' @description Given a vector of proteins, finds known proteases acting on
#' cleavage sites.
#' @param protein a vector of corresponding UniProt Accession IDs.
#' @param summarize if false (default), provides a detailed table of all
#' associated cleaving events, otherwise outputs a summarized table and only
#' includes reviewed (Uniprot) entries.
#'
#' @include Classes.R Generics.R Methods.R helper-functions.R
#'
#' @return data.table, character
#'
#' @examples
#' protein <- c("P05067", "P68871")
#' searchSubstrate(protein = protein)
#'
#' @importFrom data.table data.table
#'
#' @export

searchSubstrate <- function(protein, summarize = FALSE) {

    # Define local variables as NULL (due to non-standard evaluation
    # in data.table)
    `Protease (Uniprot)` <- `Protease (MEROPS)` <- `Substrate organism` <- NULL
    `Substrate (Uniprot)` <- `Protease status` <- .N <- NULL

    # Internal data: MEROPS Substrate_search.sql and Uniprot ID to MEROPS
    # identifier mapping
    mer <- data.table::fread(
        system.file("extdata", "mer.tab.gz", package = "proteasy"))
    merops_map <- data.table::fread(
        system.file("extdata", "merops.map.tab.gz", package = "proteasy"))

    # Show all data for substrate

    r <- mer[`Substrate (Uniprot)` %in% protein]

    r <- mapMEROPSIDs(r, merops_map)

    if(summarize == TRUE) {
        r <- r[`Protease status` == "reviewed", c(1, 7)]$`Protease (Uniprot)`
    }

    return(unique(r))

}
