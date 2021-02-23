
# needs brush up

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

packages = c("dplyr", "tidyr",
             "igraph", "Loregic",
             "reshape2","WGCNA",
             "flashClust","arsenal",
             "piano")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

require(topGO)
require(org.Hs.eg.db)
