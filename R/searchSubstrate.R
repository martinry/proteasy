
#' @name searchSubstrate
#' @title Show Cleaving Data for a Substrate by Uniprot accession
#' @usage searchSubstrate(protein, summarize = FALSE)
#' @description Given a vector of proteins, finds known proteases acting on cleavage sites.
#' @param protein a vector of corresponding UniProt Accession IDs.
#' @param summarize if false (default), provides a detailed table of all associated cleaving events, otherwise outputs a summarized table and only includes reviewed (Uniprot) entries.
#'
#' @include Classes.R Generics.R Methods.R helper-functions.R
#'
#' @return data.table, character
#'
#' @examples
#' protein <- c("P05067", "P68871")
#' searchSubstrate(protein = protein)
#'
#' @export

searchSubstrate <- function(protein, summarize = FALSE) {

    # Define local variables as NULL (due to non-standard evaluation in data.table)
    `Protease (Uniprot)` <- `Protease (MEROPS)` <- `Substrate organism` <- `Substrate (Uniprot)` <- `Protease status` <- .N <- NULL

    # Show all data for substrate

    r <- mer[`Substrate (Uniprot)` %in% protein]

    r <- mapMEROPSIDs(r)

    if(summarize) r <- r[`Protease status` == "reviewed", c(1, 7)]$`Protease (Uniprot)`

    return(unique(r))

}
