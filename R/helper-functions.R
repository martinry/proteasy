##***********************************************************************
##
##     proteasy helper functions for internal use
##
##***********************************************************************

#########################################################################
###
### Retrieve sequence data
##

getSeqData <- function(method, protein, organism) {

    if(method == "Rcpi") {

        # Rcpi method
        p <- Rcpi::getSeqFromUniProt( unique(protein), parallel = 5)

        p <- data.table::as.data.table(t(
            vapply(p, FUN.VALUE = data.table::data.table("x", "y"),
                   FUN = function(x) data.table::data.table(seq_name = as.character(
                       sub("\\|.*", "",
                           sub(".*\\|(.*)\\|.*", "\\1",
                               names(x)))), sequence = as.character(x[[1]]))
            )))

    } else if(method == "ensembldb") {

        # ensembldb method
        p <- switch(organism,
                    "Homo sapiens" = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                    "Mus musculus" = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                    "Rattus norvegicus" = EnsDb.Rnorvegicus.v79::EnsDb.Rnorvegicus.v79)

        p <- data.table::as.data.table(
            ensembldb::proteins(p,
                                filter = AnnotationFilter::UniprotFilter(protein),
                                columns = c("uniprot_id", "protein_sequence"),
                                return.type = "data.frame"))

        p <- data.table::setnames(x = p,
                                  old = c("uniprot_id", "protein_sequence"),
                                  new = c("seq_name", "sequence"))

        p <- p[!BiocGenerics::duplicated(p$seq_name), 1:2]

    }

    return(p)

}

#########################################################################
###
### Match N/C termini between user input and MEROPS db
##

matchTermini <- function(input) {

    .N <- NULL

    # Find N-terminus matches
    data.table::setkeyv(input, c("protein", "start_pos"))
    data.table::setkeyv(mer, c("Substrate (Uniprot)", "Residue number"))
    N <- input[mer, nomatch = NULL]
    N$terminus <- "N"

    # Find C-terminus matches
    data.table::setkeyv(input, c("protein", "end_pos"))
    data.table::setkeyv(mer, c("Substrate (Uniprot)", "Residue number"))
    C <- input[mer, nomatch = NULL]
    C$terminus <- "C"

    if(N[, .N] == 0 & C[, .N] == 0) stop("No matches found.", call. = FALSE)

    r <- data.table::rbindlist(list(N, C))

    return(r)


}


#########################################################################
###
### Map MEROPS ID to Uniprot
##

mapMEROPSIDs <- function(r) {

    seq_name <- .N <- `Protease organism` <- `Protease status` <- NULL

    merops_map <- merops_map[`Protease organism` %in% r$`Substrate organism`]

    data.table::setkey(r, "Protease (MEROPS)")
    data.table::setkey(merops_map, "Protease (MEROPS)")

    r <- merops_map[r, nomatch = NULL]

    r <- r[order(`Protease status`, decreasing = FALSE)]


    return(r)

}
