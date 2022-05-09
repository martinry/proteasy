
# Introduction

This package utilizes data derived from the MEROPS database (https://www.ebi.ac.uk/merops/) and facilitates retrieving proteolytic enzymes data by mapping peptide termini to known sites where a protease cleaves.

# Installation

```r
install.packages("devtools")
library("devtools")
install_github("martinry/proteasy")
```

# Usage

```r

# Example 1
searchSubstrate(protein = "P01042", summarize = TRUE)

# Example 2
searchProtease(protein = "P39900", summarize = TRUE)

# Example 3
protein <- c("P02671", "P02671", "P68871", "P01011")
peptide <- c("FEEVSGNVSPGTR", "FVSETESR", "LLVVYPW", "ITLLSAL")

res <- findProtease(protein = protein, peptide = peptide, organism = "Homo sapiens")

proteasy:::substrates(res)
proteasy:::proteases(res)
proteasy:::cleavages(res)

# Example 4

browseProtease("P07339", keytype = "UniprotID") # (opens web browser)

```

