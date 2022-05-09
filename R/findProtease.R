
#' @name findProtease
#' @title Find Proteases
#' @usage findProtease(protein, peptide, organism, start_pos, end_pos)
#' @description Given a vector of peptides and proteins, finds known proteases acting on cleavage sites.
#' @param protein a vector of UniProt Accession IDs.
#' @param peptide a vector of amino acid sequences corresponding to the proteins.
#' @param organism name of organism.
#' @param start_pos (optional) numeric vector of N-terminus positions in protein sequence.
#' @param end_pos (optional) numeric vector of C-terminus positions in protein sequence.
#'
#' @include Classes.R Generics.R Methods.R helper-functions.R
#'
#' @return S4 object Cleavages
#'
#' @examples
#' protein <- c("P02671", "P02671", "P68871", "P01011")
#' peptide <- c("FEEVSGNVSPGTR", "FVSETESR", "LLVVYPW", "ITLLSAL")
#' res <- findProtease(protein = protein, peptide = peptide, organism = "Homo sapiens")
#'
#' @importFrom data.table data.table
#'
#' @export

findProtease <- function(protein, peptide, organism = "Homo sapiens", start_pos, end_pos) {

    # Define local variables as NULL (due to non-standard evaluation in data.table)
    `Protease (Uniprot)` <- `Protease (MEROPS)` <- `Substrate organism` <- NULL
    `Protease organism` <- `Substrate (Uniprot)` <- `Protease status` <- seq_name <- .N <- NULL

    # Internal data: MEROPS Substrate_search.sql and Uniprot ID to MEROPS identifier mapping
    mer <- get0("mer", envir = asNamespace("proteasy"))
    merops_map <- get0("merops_map", envir = asNamespace("proteasy"))

    if(!(length(peptide) == length(protein))) stop("Peptide and protein vectors must be the same length.")
    if(!(organism %in% unique(mer$`Substrate organism`))) stop("Organism not recognized")

    if((missing(start_pos) | missing(end_pos))) {

        # Find start_pos and end_pos by mapping peptide against protein sequence

        if(!(organism %in% c("Homo sapiens", "Rattus norvegicus", "Mus musculus"))) {
            # Rcpi method
            p <- getSeqData(method = "Rcpi", protein = protein, organism = organism)

        } else {
            # ensembldb method
            p <- getSeqData(method = "ensembldb", protein = protein, organism = organism)
        }

        unique_proteins <- unique(protein)

        # Proteins where sequence data was not found
        unmapped <- unique_proteins[!(unique_proteins %in% p$seq_name)]

        if(length(unmapped) == length(unique_proteins)) stop("No accessions could be mapped.")

        if(length(unmapped) > 0) {
            warning(paste0(
                "Some protein accessions could not be mapped (", length(unmapped),
                "). This could be due to use of obsolete accessions or incorrect identifier type. ",
                paste0(unmapped, collapse = ", "))
                )
            }

        p <- p[seq_name %in% unique_proteins,]

        input <- unique(data.table::data.table(protein, peptide))

        data.table::setkeyv(p, "seq_name")
        data.table::setkeyv(input, "protein")

        input <- p[input]

        str_pos <- data.table::as.data.table(stringr::str_locate(input$sequence, input$peptide))

        input$start_pos <- as.character(str_pos$start)
        input$end_pos <- as.character(str_pos$end)

        names(input)[1] <- "protein"

    } else {
        # Find using position matching
        input <- data.table::data.table(protein, peptide, start_pos, end_pos)
    }

    mer <- mer[`Substrate organism` == organism & `Substrate (Uniprot)` %in% unique_proteins]

    r <- matchTermini(input)

    r <- mapMEROPSIDs(r)

    r <- r[!duplicated(r[, c("peptide", "protein", "Protease (Uniprot)",
                             "Protease (MEROPS)", "terminus", "Cleavage type")])]

    names(r) <- c("Protease (Uniprot)",     # 1
                  "Protease status",        # 2
                  "Protease organism",      # 3
                  "Protease (MEROPS)",      # 4
                  "Substrate (Uniprot)",    # 5
                  "Substrate sequence",     # 6
                  "Peptide",                # 7
                  "Start position",         # 8
                  "End position",           # 9
                  "Cleaved residue",        # 10
                  "Substrate name",         # 11
                  "Substrate organism",     # 12
                  "Protease name",          # 13
                  "Cleavage type",          # 14
                  "Cleaved terminus")       # 15

    substrate <- data.table::data.table("Substrate name"       = r[, 11][[1]],
                                        "Substrate (Uniprot)"  = r[, 5][[1]],
                                        "Substrate sequence"   = r[, 6][[1]],
                                        "Substrate length"     = nchar(r[, 6][[1]]),
                                        "Peptide"              = r[, 7][[1]],
                                        "Start position"       = r[, 8][[1]],
                                        "End position"         = r[, 9][[1]])

    protease <- data.table::data.table("Protease name"         = r[, 13][[1]],
                                       "Protease (Uniprot)"    = r[, 1][[1]],
                                       "Protease status"       = r[, 2][[1]],
                                       "Protease (MEROPS)"     = r[, 4][[1]],
                                       "Protease URL"          = paste0("https://www.ebi.ac.uk/merops/cgi-bin/pepsum?id=", r[, 4][[1]]))

    cleavage <- data.table::data.table("Substrate (Uniprot)"   = r[, 5][[1]],
                                        "Peptide"              = r[, 7][[1]],
                                        "Protease (Uniprot)"   = r[, 1][[1]],
                                        "Protease status"      = r[, 2][[1]],
                                        "Cleaved residue"      = r[, 10][[1]],
                                        "Cleaved terminus"     = r[, 15][[1]],
                                        "Cleavage type"        = r[, 14][[1]])


    return(
        methods::new(
            Class = "Cleavages",
            organism       = organism,
            substrate      = unique(substrate),
            protease       = unique(protease),
            cleavage       = cleavage
        )
    )

}

