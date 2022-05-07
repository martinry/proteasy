
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
#' @import data.table
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
            proteome <- switch(organism,
                               "Homo sapiens" = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                               "Mus musculus" = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                               "Rattus norvegicus" = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)

            proteome <-
                data.table::as.data.table(
                    ensembldb::proteins(
                        proteome,
                        filter = AnnotationFilter::UniprotFilter(protein),
                        columns = c("uniprot_id", "protein_sequence"),
                        return.type = "data.frame"
                    )
                )

            proteome <- data.table::setnames(x = proteome,
                                             old = c("uniprot_id", "protein_sequence"),
                                             new = c("seq_name", "sequence"))

            proteome <- proteome[!BiocGenerics::duplicated(proteome$seq_name), 1:2]

        }

        # List proteins not found
        unmapped <- unique_proteins[!(unique_proteins %in% proteome$seq_name)]

        if(length(unmapped) == length(unique_proteins)) stop("No accessions could be mapped.")

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

        proteome <- proteome[proteome$seq_name %in% unique_proteins,]

        input <- unique(data.table::data.table(protein, peptide))

        data.table::setkeyv(proteome, "seq_name")
        data.table::setkeyv(input, "protein")

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

    r <- data.table::rbindlist(list(N, C))

    merops_map <- merops_map[startsWith(Organism, organism)]

    data.table::setkey(r, "code")
    data.table::setkey(merops_map, "Cross-reference (MEROPS)")

    r <- merops_map[r, nomatch = NULL]

    r <- r[order(Status, decreasing = F)]

    r <- r[!duplicated(r[, c("peptide", "protein", "Entry", "Cross-reference (MEROPS)", "terminus", "cleavage_type")])]

    names(r) <- c("Protease (Uniprot)",     # 1
                  "Protease status",        # 2
                  "Protease (Gene names)",  # 3
                  "Organism",               # 4
                  "Protease (MEROPS)",      # 5
                  "Substrate (Uniprot)",    # 6
                  "Substrate sequence",     # 7
                  "Peptide",                # 8
                  "Start position",         # 9
                  "End position",           # 10
                  "Cleaved residue",        # 11
                  "Substrate name",         # 12
                  "organism",               # 13
                  "Protease name",          # 14
                  "Cleavage type",          # 15
                  "Cleaved terminus")       # 16

    substrate <- data.table::data.table("Substrate name"      = r[, 12][[1]],
                                        "Substrate (Uniprot)" = r[, 6][[1]],
                                        "Substrate sequence"  = r[, 7][[1]],
                                        "Substrate length"    = nchar(r[, 7][[1]]),
                                        "Peptide"             = r[, 8][[1]],
                                        "Start position"      = r[, 9][[1]],
                                        "End position"        = r[, 10][[1]])

    protease <- data.table::data.table("Protease name"         = r[, 14][[1]],
                                       "Protease (Uniprot)"    = r[, 1][[1]],
                                       "Protease status"       = r[, 2][[1]],
                                       "Protease (Gene names)" = r[, 3][[1]],
                                       "Protease (MEROPS)"     = r[, 5][[1]],
                                       "Protease URL"          = paste0("https://www.ebi.ac.uk/merops/cgi-bin/pepsum?id=", r[, 5][[1]]))

    cleavage <- data.table::data.table("Substrate (Uniprot)" = r[, 6][[1]],
                                        "Peptide"            = r[, 8][[1]],
                                        "Protease (Uniprot)" = r[, 1][[1]],
                                        "Protease status"    = r[, 2][[1]],
                                        "Cleaved residue"    = r[, 11][[1]],
                                        "Cleaved terminus"   = r[, 16][[1]],
                                        "Cleavage type"      = r[, 15][[1]])


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

