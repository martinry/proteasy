#' @name browseProtease
#' @title Browse Protease on MEROPS
#' @usage browseProtease(p, keytype)
#' @description Opens relevant MEROPS (https://www.ebi.ac.uk/merops/)
#' page with information on specified protease.
#' @param p a single protease
#' @param keytype UniprotID (default) or MEROPS.
#' @examples
#' \donttest{
#' if (interactive()) {
#'
#' # The following function calls opens in browser
#'
#' browseProtease("P07339", keytype = "UniprotID")
#' browseProtease("A01.009", keytype = "MEROPS")
#'
#' }
#' }
#' @return utils::browseURL
#'
#' @export

browseProtease <- function(p, keytype = "UniprotID") {

    # Define local variables as NULL (due to non-standard
    # evaluation in data.table)
    `Protease (Uniprot)` <- `Protease (MEROPS)` <- NULL

    # Internal data: MEROPS Substrate_search.sql and Uniprot ID to
    # MEROPS identifier mapping
    merops_map <- data.table::fread(
        system.file("extdata", "merops.map.tab.gz", package = "proteasy"))

    if(keytype == "UniprotID"){
        if(!(p %in% merops_map$`Protease (Uniprot)`)) {
            stop("Identified not recognized.")
            }

        p <- merops_map[`Protease (Uniprot)` == p]$`Protease (MEROPS)`
    } else if(keytype == "MEROPS") {
        if(!(p %in% merops_map$`Protease (MEROPS)`)) {
            stop("Identified not recognized.")
        } 
    } else {
        stop('Not a valid keytype. Valid keytypes are "UniprotID", "MEROPS".')
    }

    utils::browseURL(
        paste0("https://www.ebi.ac.uk/merops/cgi-bin/pepsum?id=", p)
        )

}
