
#' @name findProtease
#' @title Find Proteases
#' @usage findProtease(peptide, protein)
#' @description Given a vector of peptides and proteins, finds proteases acting on cleavage sites
#' @param protein a vector of corresponding UniProt Accession IDs
#' @param peptide a vector of amino acid sequences
#' @param organism  name of organism
#' @param start_pos  (optional) vector of N-terminus position in protein sequence
#' @param end_pos  (optional) vector of C-terminus position in protein sequence
#'
#' @return S4 object proteasy_res
#'
#' @examples
#' protein <- c("P02671", "P02671", "P68871", "P01011")
#' peptide <- c("FEEVSGNVSPGTR", "FVSETESR", "LLVVYPW", "ITLLSAL")
#' res <- findProtease(protein = protein, peptide = peptide)
#'
#'
#' @export

findProtease <- function(protein, peptide,
                         start_pos, end_pos,
                         organism = "Homo sapiens") {

    if(!(length(peptide) == length(protein))) stop("Peptide and protein vectors must be the same length.")
    if(!(organism %in% unique(mer$organism))) stop("Organism not recognized")

    if((missing(start_pos) | missing(end_pos))) {

        unique_proteins <- unique(protein)

        # List proteins not found
        unmapped <- unique_proteins[!(unique_proteins %in% proteome$seq_name)]

        if(length(unmapped) > 0) {
            warning(
                paste0(
                    "Some protein accessions could not be mapped (",
                    length(unmapped),
                    "). This could be due to use of obsolete accessions or incorrect identifier type. ",
                    paste0(unmapped, collapse = ", ")
                )
            )
        }

        proteome <- proteome[seq_name %in% unique_proteins]

        input <- unique(data.table(protein, peptide))

        data.table::setkey(proteome, "seq_name")
        data.table::setkey(input, "protein")

        input <- proteome[input]

        str_pos <- stringr::str_locate(input$sequence, input$peptide) %>% as.data.table

        input$start_pos <- as.character(str_pos$start)
        input$end_pos <- as.character(str_pos$end)

        names(input)[1] <- "protein"

    } else {
        # Find using position matching
        input <- data.table::data.table(protein, peptide, start_pos, end_pos)
    }

    mer <- mer[organism == organism]
    mer <- mer[Uniprot %in% unique_proteins]

    # Find N-terminus matches
    data.table::setkeyv(input, c("protein", "start_pos"))
    data.table::setkeyv(mer, c("Uniprot", "resnum"))
    N <- input[mer, nomatch = NULL]
    N$terminus <- "N"

    # Find C-terminus matches
    data.table::setkeyv(input, c("protein", "end_pos"))
    data.table::setkeyv(mer, c("Uniprot", "resnum"))
    C <- input[mer, nomatch = NULL]
    C$terminus <- "C"

    if(N[, .N] == 0 & C[, .N] == 0) stop("No matches found.")

    res <- rbindlist(list(N, C))

    merops_map <- merops_map[startsWith(Organism, organism)]

    setkey(res, "code")
    setkey(merops_map, "Cross-reference (MEROPS)")

    res2 <- merops_map[res, nomatch = NULL]

    res2 <- res2[order(Status, decreasing = F)]

    res2 <- res2[!duplicated(res2[, c("peptide", "protein", "Entry", "Cross-reference (MEROPS)", "terminus", "cleavage_type")])]

    names(res2) <-
        c(
            "Protease (Uniprot)",
            "Protease status",
            "Protease (Gene names)",
            "Organism",
            "Protease (MEROPS)",
            "Substrate (Uniprot)",
            "Substrate sequence",
            "Peptide",
            "Start position",
            "End position",
            "Cleaved residue",
            "Substrate name",
            "organism",
            "Protease name",
            "Cleavage type",
            "Cleaved terminus"
        )

    substrate <- data.table("Substrate name" = res2$`Substrate name`,
                            "Substrate (Uniprot)" = res2$`Substrate (Uniprot)`,
                            "Substrate sequence" = res2$`Substrate sequence`,
                            "Substrate length" = nchar(res2$`Substrate sequence`),
                            "Peptide" = res2$Peptide,
                            "Start position" = res2$`Start position`,
                            "End position" = res2$`End position`)

    protease <- data.table("Protease name" = res2$`Protease name`,
                           "Protease (Uniprot)" = res2$`Protease (Uniprot)`,
                           "Protease status" = res2$`Protease status`,
                           "Protease (Gene names)" = res2$`Protease (Gene names)`,
                           "Protease (MEROPS)" = res2$`Protease (MEROPS)`,
                           "Protease URL" = paste0("https://www.ebi.ac.uk/merops/cgi-bin/pepsum?id=", res2$`Protease (MEROPS)`))

    output <- data.table("Substrate (Uniprot)" = res2$`Substrate (Uniprot)`,
                         "Peptide" = res2$Peptide,
                         "Protease (Uniprot)" = res2$`Protease (Uniprot)`,
                         "Cleaved residue" = res2$`Cleaved residue`,
                         "Cleaved terminus" = res2$`Cleaved terminus`,
                         "Cleavage type" = res2$`Cleavage type`
                         )


    res <- new("proteasy_res",
               organism       = organism,
               substrate      = unique(substrate),
               protease       = unique(protease),
               output         = output
    )

    return(res)
}

