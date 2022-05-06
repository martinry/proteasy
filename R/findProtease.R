
#' @name findProtease
#' @title Find Proteases
#' @usage findProtease(protein, peptide, organism, start_pos, end_pos)
#' @description Given a vector of peptides and proteins, finds proteases acting on cleavage sites
#' @param protein a vector of corresponding UniProt Accession IDs
#' @param peptide a vector of amino acid sequences
#' @param organism  name of organism
#' @param start_pos  (optional) vector of N-terminus position in protein sequence
#' @param end_pos  (optional) vector of C-terminus position in protein sequence
#'
#' @include AllClasses.R methods.R helper-functions.R
#'
#' @return S4 object Cleavages
#'
#' @examples
#' protein <- c("P02671", "P02671", "P68871", "P01011")
#' peptide <- c("FEEVSGNVSPGTR", "FVSETESR", "LLVVYPW", "ITLLSAL")
#' res <- findProtease(protein = protein, peptide = peptide, organism = "Homo sapiens")
#'
#' @export

findProtease <- function(protein, peptide, organism = "Homo sapiens", start_pos, end_pos) {

    # Define local variables as NULL (due to non-standard evaluation in data.table)
    uniprot_id <- seq_name <- Uniprot <- Organism <- Status <- .N <- NULL

    mer <- get0("mer", envir = asNamespace("proteasy"))
    merops_map <- get0("merops_map", envir = asNamespace("proteasy"))

    if(!(length(peptide) == length(protein))) stop("Peptide and protein vectors must be the same length.")
    if(!(organism %in% unique(mer$organism))) stop("Organism not recognized")

    if((missing(start_pos) | missing(end_pos))) {

        unique_proteins <- unique(protein)

        if(!(organism %in% c("Homo sapiens", "Rattus norvegicus", "Mus musculus"))) {

            # Retrieve sequences using Rcpi
            proteome <- Rcpi::getSeqFromUniProt( unique_proteins, parallel = 5)

            proteome <-
                data.table::as.data.table(t(mapply(
                    proteome,
                    FUN = function(x)
                        data.table::data.table(seq_name = as.character(sub(
                            "\\|.*", "", sub(".*\\|(.*)\\|.*", "\\1", names(x)))
                        ), sequence = as.character(x[[1]]))
                )))

        } else {
            proteome <- orgDB(organism)

            proteome <-
                data.table::as.data.table(
                    ensembldb::proteins(
                        proteome,
                        filter = AnnotationFilter::UniprotFilter(protein),
                        columns = c("uniprot_id", "protein_sequence"),
                        return.type = "data.frame"
                    )
                )

            proteome <- proteome[!duplicated(uniprot_id), 1:2]
            names(proteome) <- c("seq_name", "sequence")

        }

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

        input <- unique(data.table::data.table(protein, peptide))

        data.table::setkey(proteome, "seq_name")
        data.table::setkey(input, "protein")

        input <- proteome[input]

        str_pos <- data.table::as.data.table(stringr::str_locate(input$sequence, input$peptide))

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

    res <- data.table::rbindlist(list(N, C))

    merops_map <- merops_map[startsWith(Organism, organism)]

    data.table::setkey(res, "code")
    data.table::setkey(merops_map, "Cross-reference (MEROPS)")

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

    substrate <- data.table::data.table("Substrate name" = res2$`Substrate name`,
                            "Substrate (Uniprot)" = res2$`Substrate (Uniprot)`,
                            "Substrate sequence" = res2$`Substrate sequence`,
                            "Substrate length" = nchar(res2$`Substrate sequence`),
                            "Peptide" = res2$Peptide,
                            "Start position" = res2$`Start position`,
                            "End position" = res2$`End position`)

    protease <- data.table::data.table("Protease name" = res2$`Protease name`,
                           "Protease (Uniprot)" = res2$`Protease (Uniprot)`,
                           "Protease status" = res2$`Protease status`,
                           "Protease (Gene names)" = res2$`Protease (Gene names)`,
                           "Protease (MEROPS)" = res2$`Protease (MEROPS)`,
                           "Protease URL" = paste0("https://www.ebi.ac.uk/merops/cgi-bin/pepsum?id=", res2$`Protease (MEROPS)`))

    cleavage <- data.table::data.table("Substrate (Uniprot)" = res2$`Substrate (Uniprot)`,
                         "Peptide" = res2$Peptide,
                         "Protease (Uniprot)" = res2$`Protease (Uniprot)`,
                         "Protease status" = res2$`Protease status`,
                         "Cleaved residue" = res2$`Cleaved residue`,
                         "Cleaved terminus" = res2$`Cleaved terminus`,
                         "Cleavage type" = res2$`Cleavage type`
                         )


    res <- methods::new(
        Class = "Cleavages",
        organism       = organism,
        substrate      = unique(substrate),
        protease       = unique(protease),
        cleavage       = cleavage
    )

    return(res)
}

