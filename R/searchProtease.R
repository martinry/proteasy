
#' @name searchProtease
#' @title Show Cleaving Data for a Peptidase or Inhibitor by Uniprot accession
#' @usage searchProtease(protein, summarize = FALSE)
#' @description Given a vector of proteins, finds which substrates they cleave.
#' @param protein a vector of corresponding UniProt Accession IDs.
#' @param summarize if false (default), provides a detailed table of all associated cleaving events, otherwise outputs a summarized table and only includes reviewed (Uniprot) entries.
#'
#' @include Classes.R Generics.R Methods.R helper-functions.R
#'
#' @return data.table, character
#'
#' @examples
#' protein <- c("P98073", "P00734")
#' searchProtease(protein = protein)
#'
#' @export
#'
#' @importFrom data.table data.table

searchProtease <- function(protein, summarize = FALSE) {

    # Define local variables as NULL (due to non-standard evaluation in data.table)
    `Protease (Uniprot)` <- `Protease (MEROPS)` <- `Substrate organism` <- `Substrate (Uniprot)` <- `Protease status` <- .N <- NULL

    # Show all data for protease

    protein_merops <- merops_map[`Protease (Uniprot)` %in% protein, 3:4]

    r <- mer[`Protease (MEROPS)` %in% protein_merops$`Protease (MEROPS)` & `Substrate (Uniprot)` != "\\N" & `Substrate organism` %in% protein_merops$`Protease organism`]

    r <- mapMEROPSIDs(r)

    if(summarize) r <- r[`Protease status` == "reviewed"]$`Substrate (Uniprot)`

    return(unique(r))

}
