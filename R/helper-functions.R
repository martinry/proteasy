orgDB <- function(s) {
    switch(s,
           "Homo sapiens" = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
           "Mus musculus" = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
           "Rattus norvegicus" = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
}
