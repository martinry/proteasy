
#' @name searchProtease
#' @title Show Cleaving Data for a Peptidase or Inhibitor by Uniprot accession
#' @usage searchProtease(protein, organism = "Homo sapiens", summarize = FALSE)
#' @description Given a vector of proteins, finds which substrates they cleave.
#' @param protein a vector of corresponding UniProt Accession IDs.
#' @param organism name of organism.
#' @param summarize if false (default), provides a detailed table of all
#' associated cleaving events, otherwise outputs a summarized table and only
#' includes reviewed (Uniprot) entries.
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
#' @import R.utils
#'
#' @export

searchProtease <- function(protein,
                           organism = "Homo sapiens",
                           summarize = FALSE) {

    # Define local variables as NULL (due to non-standard
    # evaluation in data.table)
    `Protease (Uniprot)` <- `Protease (MEROPS)` <- NULL
    `Substrate (Uniprot)` <- `Protease status`  <- .N <- NULL
    `Protease organism` <- `Substrate organism` <- NULL

    # Internal data: MEROPS Substrate_search.sql and
    # Uniprot ID to MEROPS identifier mapping
    mer <- data.table::fread(
        system.file("extdata", "mer.tab.gz", package = "proteasy"))
    merops_map <- data.table::fread(
        system.file("extdata", "merops.map.tab.gz", package = "proteasy"))

    # Show all data for protease

    protein_merops <- merops_map[`Protease (Uniprot)` %in% protein &
                                     `Protease organism` == organism]
    
    data.table::setkeyv(protein_merops, "Protease (MEROPS)")
    data.table::setkeyv(mer, "Protease (MEROPS)")
    
    r <- mer[protein_merops, nomatch = NULL]
    
    r <- r[`Substrate organism` == organism]

    if(summarize == TRUE) {
        r <- r[`Protease status` == "reviewed"]$`Substrate (Uniprot)`
    }

    return(unique(r))

}
