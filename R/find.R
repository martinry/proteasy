
#' @name find
#' @usage find(peptide, protein)
#' @description This is an allometric function to return the tree volume
#' @param peptide a vector of amino acid sequences
#' @param protein a vector of corresponding UniProt Accession IDs
#' @param start_pos  (optional) vector of N-terminus position in protein sequence
#' @param end_pos  (optional) vector of C-terminus position in protein sequence
#'
#' @return matching_cleavages matching_cleavages
#' @export

find <- function(peptide,
                 protein,
                 start_pos,
                 end_pos,
                 organism = "Homo sapiens") {

    unique_proteins <- protein %>% unique

    mer <- fread("../lib/Substrate_search.txt")

    mer <- mer[organism == organism]

    if((missing(start_pos) | missing(end_pos))) {

        # Find by aligning w/ proteome

        library("Biostrings", quietly = T)

        proteome <- Biostrings::readAAStringSet("../lib/uniprot-proteome_UP000005640+reviewed_yes.fasta.gz")
        seq_name = names(proteome)
        sequence = paste(proteome)
        df <- data.table(seq_name, sequence)
        df$seq_name %<>% sub(".*\\|(.*)\\|.*", "\\1", .)
        df$seq_name %<>% sub("\\|.*", "\\1", .)

        # List proteins not found
        unmapped <- unique_proteins[!(unique_proteins %in% df$seq_name)]

        if(length(unmapped) > 0) warning(paste0("Some protein accessions could not be mapped (", length(unmapped),
                                                "). This could be due to use of obsolete accessions or incorrect identifier type. ",
                                                unmapped))

        df <- df[seq_name %in% unique_proteins]

        input <- data.table(protein, peptide)

        setkey(df, "seq_name")
        setkey(input, "protein")

        input <- df[input]

        str_pos <- str_locate(input$sequence, input$peptide) %>% as.data.table

        input$start_pos <- str_pos$start %>% as.character
        input$end_pos <- str_pos$end %>% as.character

        names(input)[1] <- "protein"



    } else {

        # Find using position matching

        input <- data.table(protein, peptide, start_pos, end_pos)


    }

    mer <- mer[Uniprot %in% unique_proteins]
    mer$resnum %<>% as.character

    matching_cleavages <- data.table()


    for(i in unique_proteins) {

        mp <- mer[Uniprot == i]
        lfqp <- input[protein == i]

        mp <- mp[resnum %in% lfqp$start | resnum %in% lfqp$end]

        start <- mp[resnum %in% lfqp$start]
        end <- mp[resnum %in% lfqp$end]

        if(start[, .N] > 0) {
            start_cleavages <- data.table("Substrate" = start$Uniprot,
                                          "ProteaseName" = start$Protease,
                                          "Residue" = start$Letter,
                                          "Resnum" = start$resnum,
                                          "Terminus" = "N")

            setkey(start_cleavages, "Resnum")
            setkey(lfqp, "start_pos")

            start_cleavages <- start_cleavages[lfqp, nomatch = 0 , allow.cartesian = T]
            start_cleavages <- start_cleavages[, c("Substrate", "ProteaseName", "Residue", "Resnum", "Terminus", "peptide")]


        }


        if(end[, .N] > 0) {
            end_cleavages <- data.table("Substrate" = end$Uniprot,
                                        "ProteaseName" = end$Protease,
                                        "Residue" = end$Letter,
                                        "Resnum" = end$resnum,
                                        "Terminus" = "C")

            setkey(end_cleavages, "Resnum")
            setkey(lfqp, "end_pos")

            end_cleavages <- end_cleavages[lfqp, nomatch = 0, allow.cartesian = T]
            end_cleavages <- end_cleavages[, c("Substrate", "ProteaseName", "Residue", "Resnum", "Terminus", "peptide")]
        }

        if(mp[, .N] > 0) {
            cleavages <- rbindlist(list(start_cleavages, end_cleavages))

            matching_cleavages <- rbindlist(list(matching_cleavages, cleavages))
        }




    }

    return(matching_cleavages)



}
