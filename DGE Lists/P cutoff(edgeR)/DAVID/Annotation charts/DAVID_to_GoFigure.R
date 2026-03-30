library(dplyr)
filenames <- list.files(".", pattern = ".txt")
DAVID.tables <- lapply(filenames, read.table, 
                       header = TRUE, sep = "\t", quote = "")
names(DAVID.tables) <- sub(".txt", "", filenames)
names(DAVID.tables) <- gsub(" ", "", names(DAVID.tables))
for (i in 1:length(DAVID.tables)) {
  DAVID.tables[[i]]$Term <- gsub("\"", "", DAVID.tables[[i]]$Term)
  }
row_num <- length(DAVID.tables)
DAVID.KEGG <- vector(mode = "list", length = row_num)
DAVID.KEGG <- lapply(DAVID.tables, filter, Category=="KEGG_PATHWAY")

DAVID.GO <- vector(mode = "list", length = row_num)
DAVID.GO <- lapply(DAVID.tables, filter, Category!="KEGG_PATHWAY")
DAVID.GO <- lapply(DAVID.GO, dplyr::select, Term, PValue, Count)

for (i in 1:length(DAVID.GO)){
  DAVID.GO[[i]]$Term <- unlist(lapply(DAVID.GO[[i]]$Term, substr, 1, 10))
  write.table(DAVID.GO[[i]], file = paste("./GO-Figure/", names(DAVID.GO)[i], ".tsv" 
                                          ,sep = ""),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}
remove(filenames, row_num, DAVID.tables)

