library(S4Vectors)
library(edgeR)
library(dplyr)
library("AnnotationHub")
library(DESeq2)
library(gplots)
library(ggplot2)
library(clusterProfiler) #Warning: Masks select() from AnnotationDbi

############## Initial read-in #################################################
#Note: UWBC gave us an RSEM output file with the first column starting with
#r109-M, rather than a header for the gene symbol column. The entire first row
#needed to be manually shifted over by 1 to align the headers with their data,
#and "Gene Symbol" was inserted.

# Read in sample table
samples <- read.delim("SampleTable.txt", stringsAsFactors = TRUE)
samples <-samples[grep("128.M", samples$Sample, invert = TRUE),] #128M is outlier by large margin
samples$Sample <- unfactor(samples$Sample)
samples$group <- paste(samples$treatment, 
                       samples$sex, 
                       samples$location, 
                       sep="_") |> factor()
groups <- paste(samples$group) |> factor()

# Read in count table
GenewiseCounts <- read.delim("rsem_GENE(shifted).txt", row.names = "Gene.Symbol")
GenewiseCounts <- GenewiseCounts[,grep("128.M", colnames(GenewiseCounts), invert = TRUE)]
GenewiseCounts <- GenewiseCounts[,samples$Sample]


######################### edgeR ################################################
############## Load into edgeR #################################################

zb <- DGEList(GenewiseCounts, group = groups, genes = rownames(GenewiseCounts))
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

############## Filter and normalize#############################################
#Filtering out low-count genes: a common threshold is [ > 10/(minimum library
#size in millions)] in at least [minimum number of replicates] libraries. Given
#we have 48 libraries, this is a very generous threshold.
#Another approach is to use edgeR's filterByExpr() method. It has a similar
#approach, but has a more involved statistical model to determine cutoffs.

zbkeep <- filterByExpr(zb, design)
zbnorm <- zb[zbkeep, , keep.lib.sizes = FALSE]
zbnorm <- calcNormFactors(zbnorm, method = "TMM")

paste("Initial genes:", nrow(zb$genes),"Remaining genes:", nrow(zbnorm$genes))

# plot(zbnorm$samples$norm.factors, main = "Normalization factors")

#Dispersion estimation
zbnorm <- estimateDisp(zbnorm, design, robust=TRUE)

# plotBCV(zbnorm, yaxp = c(0, 3, 6))

############## Multidimensional scaling (MDS) plot #############################
colors <- c()
for (i in 1:length(samples$Sample)) {
  comp <- samples$treatment[i]
  if (comp == "Paired"){colors <- c(colors, "darkgreen")
  } else if (comp == "Separated") {colors <- c(colors, "cyan3")
  } else if (comp == "Unpaired") {colors <- c(colors, "darkorange")
  } else {print( paste0( samples$Sample[i], " is not grouped correctly!"))}
}

pch <- c()
for (i in 1:length(samples$Sample)) {
  comp <- paste(samples$sex[i], samples$location[i], sep = ", ")
  if (comp == "Male, Caudal"){pch <- c(pch, "15")
  } else if (comp == "Male, Rostral") {pch <- c(pch, "18")
  } else if (comp == "Female, Caudal") {pch <- c(pch, "16")
  } else if (comp == "Female, Rostral") {pch <- c(pch, "17")}
  else {print(paste0(samples$Sample[i], "is not grouped correctly!"))}
}
pch <- as.numeric(pch)
plotMDS(zbnorm, pch = pch, col = colors, main = "MDS plot, top 500")
legend("center",
       pch = 19,
       legend = levels(factor(samples$treatment)),
       col = c("darkgreen", "cyan3","darkorange"),
       title = "Treatment")
remove(colors, comp, pch)

############## PCA plot ########################################################
# With ggbiplot and prcomp
pca_analysis.cpm <- prcomp(t(cpm(zbnorm, log=TRUE))) #This works via cpm...Alternatives?
summary(pca_analysis.cpm)

# All groups delineated
ggbiplot.PCA <- ggbiplot::ggbiplot(pca_analysis.cpm,
                                   groups = paste(samples$treatment, samples$sex, samples$location, sep=", ") |>
                                     factor(), 
                                   ellipse = FALSE,
                                   var.axes = FALSE,
                                   scale = 0,
                                   labels = NULL,
                                   choices = c(1,2)) #Change numbers to plot other PCs
ggbiplot.PCA + theme_bw() + 
  theme(legend.title = element_blank(),
        legend.background = element_blank())

# Sex and location
ggbiplot.PCA <- ggbiplot::ggbiplot(pca_analysis.cpm,
                                   groups = paste(samples$sex, samples$location, sep=", ") |>
                                     factor(), 
                                   ellipse = FALSE,
                                   var.axes = FALSE,
                                   scale = 0,
                                   ellipse.fill = FALSE,
                                   labels = NULL,
                                   choices = c(1,2)) #Change numbers to plot other PCs
ggbiplot.PCA + theme_bw() + 
  theme(legend.position = "inside", 
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black")
  )
remove(pca_analysis.cpm, ggbiplot.PCA)
############## Generate MD plots ###############################################
for(i in 1:ncol(zbnorm)){
  plotMD(zbnorm, column=i)
  abline(h=0, col="red", lty=2, lwd=2)
}
############## Contrasts #######################################################
contrasts <- makeContrasts(
  #Comparison by treatment ACROSS groups
  "Caudal Female, Paired vs Unpaired"  = Paired_Female_Caudal  - Unpaired_Female_Caudal,            
  "Rostral Female, Paired vs Unpaired" = Paired_Female_Rostral - Unpaired_Female_Rostral,           
  "Caudal Female, Separated vs Unpaired"  = Separated_Female_Caudal  - Unpaired_Female_Caudal,         
  "Rostral Female, Separated vs Unpaired" = Separated_Female_Rostral - Unpaired_Female_Rostral,         
  "Caudal Female, Paired vs Separated"  = Paired_Female_Caudal  - Separated_Female_Caudal, 
  "Rostral Female, Paired vs Separated" = Paired_Female_Rostral - Separated_Female_Rostral,
  
  "Caudal Male, Paired vs Unpaired"  = Paired_Male_Caudal  - Unpaired_Male_Caudal,               
  "Rostral Male, Paired vs Unpaired" = Paired_Male_Rostral - Unpaired_Male_Rostral,              
  "Caudal Male, Separated vs Unpaired"  = Separated_Male_Caudal  - Unpaired_Male_Caudal,            
  "Rostral Male, Separated vs Unpaired" = Separated_Male_Rostral - Unpaired_Male_Rostral,           
  "Caudal Male, Paired vs Separated"  = Paired_Male_Caudal  - Separated_Male_Caudal,    
  "Rostral Male, Paired vs Separated" = Paired_Male_Rostral - Separated_Male_Rostral,   
  
  #Comparison by region WITHIN groups
  "Paired Female, Rostral vs Caudal"  = Paired_Female_Rostral  - Paired_Female_Caudal,   
  "Paired Male, Rostral vs Caudal"  = Paired_Male_Rostral  - Paired_Male_Caudal,   
  "Separated Female, Rostral vs Caudal"  = Separated_Female_Rostral  - Separated_Female_Caudal,
  "Separated Male, Rostral vs Caudal"  = Separated_Male_Rostral  - Separated_Male_Caudal,
  "Unpaired Female, Rostral vs Caudal"  = Unpaired_Female_Rostral  - Unpaired_Female_Caudal, 
  "Unpaired Male, Rostral vs Caudal"  = Unpaired_Male_Rostral  - Unpaired_Male_Caudal, 
  
  #Pooled comparisons
  "Pooled Male, Rostral vs Caudal" = (Paired_Male_Rostral + Separated_Male_Rostral + Unpaired_Male_Rostral) - 
    (Paired_Male_Caudal + Separated_Male_Caudal + Unpaired_Male_Caudal),     
  "Pooled Female, Rostral vs Caudal" = (Paired_Female_Rostral + Separated_Female_Rostral + Unpaired_Female_Rostral) - 
    (Paired_Female_Caudal + Separated_Female_Caudal + Unpaired_Female_Caudal),     
  "Pooled Male and Female, Rostral vs Caudal" = (Paired_Male_Rostral + Separated_Male_Rostral + Unpaired_Male_Rostral) + 
    (Paired_Female_Rostral + Separated_Female_Rostral + Unpaired_Female_Rostral) - 
    (Paired_Male_Caudal + Separated_Male_Caudal + Unpaired_Male_Caudal) - 
    (Paired_Female_Caudal + Separated_Female_Caudal + Unpaired_Female_Caudal),     
  "Pooled Caudal, Female vs Male"  = (Paired_Female_Caudal + Separated_Female_Caudal + Unpaired_Female_Caudal) - 
    (Paired_Male_Caudal + Separated_Male_Caudal + Unpaired_Male_Caudal),     
  "Pooled Rostral, Female vs Male" = (Paired_Female_Rostral + Separated_Female_Rostral + Unpaired_Female_Rostral) - 
    (Paired_Male_Rostral + Separated_Male_Rostral + Unpaired_Male_Rostral),  
  "Pooled Rostral and Caudal, Female vs Male" = (Paired_Female_Caudal + Separated_Female_Caudal + Unpaired_Female_Caudal) + 
    (Paired_Female_Rostral + Separated_Female_Rostral + Unpaired_Female_Rostral) - 
    (Paired_Male_Caudal + Separated_Male_Caudal + Unpaired_Male_Caudal) - 
    (Paired_Male_Rostral + Separated_Male_Rostral + Unpaired_Male_Rostral),  
  
  levels=design
)
############## edgeR DGE list generating and contrast MD plots #################
fit <- glmQLFit(zbnorm, design)
plotQLDisp(fit)

#Heat map setup
logCPM <- cpm(zbnorm, log = TRUE, normalized.lib.sizes = TRUE)
rownames(logCPM) <- rownames(zbnorm$genes)
colnames(logCPM) <- paste(zbnorm$samples$group, 1:2, sep = "-")

#DGE analysis loop setup
edgeR_list <- list() # This will become a list of topTags results. Each column is a contrast.
edgeR_list_pcutoff <- list()
edgeR_list_GO_up <- list()
edgeR_list_GO_down <- list()

for(comparison in 1:ncol(contrasts)){
  con <- colnames(contrasts)[comparison]
  tr <- glmQLFTest(fit, contrast = contrasts[, comparison])
  print(paste(comparison, ":",colnames(contrasts)[comparison]))
  
  # DE plots
  png(filename = paste("./DEG Plots(edgeR)/",
                       con, ".png",
                       sep = ""),
      height = 400, width = 400)
  plotMD(tr, legend = "topright", main = paste(con))
  dev.off()
  
  # Remainder of loop relates to writing DEG lists and related files
  if (comparison == 1){is_de <- decideTests(tr)
  } else {is_de <- cbind(is_de, decideTests(tr))}
  
  edgeR_list[[comparison]] <- topTags(tr, n = Inf)[[1]]
  names(edgeR_list)[comparison] <-  con
  
  # List of q <0.05
  edgeR_list_pcutoff[[comparison]] <- edgeR_list[[comparison]][edgeR_list[[comparison]]$FDR < 0.05, ]
  names(edgeR_list_pcutoff)[comparison] <-  con
  
  # Separated by sign and LFC threshold = 1 (Tested 0, 1, 1.5 and 2 on Pooled R-C. Fewest unique/abberant GO categories in 1.5)
  edgeR_list_GO_up[[comparison]] <- edgeR_list[[comparison]][(edgeR_list[[comparison]]$FDR < 0.05 & edgeR_list[[comparison]]$logFC > 1.5),]
  names(edgeR_list_GO_up)[comparison] <-  con
  
  edgeR_list_GO_down[[comparison]] <- edgeR_list[[comparison]][(edgeR_list[[comparison]]$FDR < 0.05 & edgeR_list[[comparison]]$logFC < -1.5),]
  names(edgeR_list_GO_down)[comparison] <-  con
  
  #Heat map
  o <- order(tr$table$PValue) #Heat map setup is contrast-dependent
  logCPM_sub <- logCPM[o[1:100], ]
  logCPM_sub <- t(scale(t(logCPM_sub)))
  col.pan <- colorpanel(100, "blue", "white", "red")
  # png(filename = paste("./Heatmaps(edgeR)/",
  #                      con, ".png",
  #                      sep = ""),
  #     height = 800, width = 1200)
  # heatmap.2(logCPM_sub,
  #           col = col.pan,
  #           Rowv = TRUE,
  #           scale = "none",
  #           trace = "none",
  #           dendrogram = "both",
  #           density.info = "none",
  #           key = FALSE,
  #           margins = c(8,8),
  #           main = colnames(contrasts)[comparison])
  # dev.off()
}
colnames(is_de) <- colnames(contrasts)
# Add gene names to data frame:
tgu <- query(AnnotationHub(), c("Taeniopygia guttata", "OrgDb"))[["AH115578"]]
for(comparison in 1:ncol(contrasts)){
  edgeR_list_pcutoff[[comparison]][7] <- AnnotationDbi::select(tgu,
                                                               keys = rownames(edgeR_list_pcutoff[[comparison]]),
                                                               columns = "CHR",
                                                               keytype = "SYMBOL")[2]
  edgeR_list_pcutoff[[comparison]][8] <- AnnotationDbi::select(tgu,
                                                               keys = rownames(edgeR_list_pcutoff[[comparison]]),
                                                               columns = "GENENAME",
                                                               keytype = "SYMBOL")[2]
}
# Lists to CSV
for(comparison in 1:ncol(contrasts)){
  con <- colnames(contrasts)[comparison]
  #Full CSV
  write.csv(edgeR_list_pcutoff[[comparison]],
            paste("./DGE Lists/P cutoff(edgeR)/",
                  con,
                  ".csv",
                  sep = "")
  )
  write.csv(edgeR_list[[comparison]],
            paste("./DGE Lists/No cutoff(edgeR)/",
                  con,
                  ".csv",
                  sep = "")
  )
  #Gene names, upregulated
  write.table(edgeR_list_GO_up[[comparison]][,1],
              paste("./DGE Lists/P cutoff(edgeR)/DAVID/",
                    con,
                    "-UP.txt",
                    sep = ""),
              sep = "\t",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE
  )
  #Gene names, downregulated
  write.table(edgeR_list_GO_down[[comparison]][,1],
              paste("./DGE Lists/P cutoff(edgeR)/DAVID/",
                    con,
                    "-DOWN.txt",
                    sep = ""),
              sep = "\t",
              col.names = FALSE,
              row.names = FALSE, 
              quote = FALSE
  )
  #Gene names, up- and downregulated
  write.table(c(edgeR_list_GO_up[[comparison]][,1], edgeR_list_GO_down[[comparison]][,1]),
              paste("./DGE Lists/P cutoff(edgeR)/DAVID/",
                    con,
                    "-ALL.txt",
                    sep = ""),
              sep = "\t",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE
  )
}
write.table(edgeR_list[[1]][,1],
            "./DGE Lists/P cutoff(edgeR)/DAVID/DAVIDBackgroundGeneList.txt",
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)
# p-value plots, for quality control
for(comparison in 1:ncol(contrasts)){
  hist(edgeR_list[[comparison]][,5], 
       main = paste(colnames(contrasts)[comparison], ", edgeR", sep = ""), 
       xlab = "p-value")
}
remove(col.pan, comparison, con, tr)

############### Venn diagrams ######################
#This section is very bloated- it is functioning as a dumping ground for figures
# for presentations/etc.

# Template
vennDiagram(is_de[,c(23,22,24)], #accepts up to 5 contrasts 
            include = c("up","down"),
            circle.col = c("red","blue","orange","green3"),
            counts.col = c("red","blue")
)
# Treatment- 4 way
# Females
vennDiagram(is_de[,c(7,8,9,10)],
            include = c("up","down"),
            circle.col = c("red","blue","orange","green3"),
            counts.col = c("red","blue")
)
# Males
vennDiagram(is_de[,c(13,14,15,16)],
            include = c("up","down"),
            circle.col = c("red","blue","orange","green3"),
            counts.col = c("red","blue")
)
# Sex differences
vennDiagram(is_de[,c(22, 23)],
            include = c("up","down"),
            circle.col = c("orange","green4"),
            counts.col = c("hotpink","blue"),
            names = c("Caudal", "Rostral"),
            show.include = "FALSE"
)
# Rostrocaudal comparison
vennDiagram(is_de[,c(20, 19)],
            include = c("up","down"),
            circle.col = c("hotpink","blue"),
            counts.col = c("green4","orange"),
            names = c("Females", "Males"),
            show.include = "FALSE"
)

############## GO enrichment and GSEA ##########################################
# There is an http issue with KEGG that stems from old versions of clusterProfiler.
# Download from GitHub if this pops up.
# remotes::install_github("YuLab-SMU/clusterProfiler") 
#

#Optional: check for up-to-date annotation package
#Search for all T. guttata DBs:
# tgutDBs <- query(AnnotationHub(), c("Taeniopygia guttata","OrgDb"))
# print(tgutDBs)
# tgu <- tgutDBs[["AH115578"]] #The specified DB was the only one found. Most recent check: 11/11/24

#Direct (ie, skipping search)
gene_background <- select(tgu,
                          keys = unlist(zbnorm$genes),
                          columns = "GID", 
                          keytype = "SYMBOL")

gene_background$SYMBOL[duplicated(gene_background$SYMBOL)]
# Four gene symbols correspond to more than one GID in NCBI. The duplicates have
# to be removed manually based on version information on NCBI (there is no consistency
# to whether the up-to-date version is first). 
gene_background <- gene_background[-c(869,1877,2217,3854),]

# Verify duplicates are gone: 
length(gene_background$SYMBOL[duplicated(gene_background$SYMBOL)]) == 0 #should be TRUE
nrow(gene_background) == nrow(zbnorm$genes) #Should be TRUE

GID_background <- unlist(gene_background[!is.na(gene_background$GID),][2])
gene_background <- unlist(gene_background[1])

DE_list_up <- list()
DE_list_down <- list()

for(i in 1:ncol(contrasts)){
  DE_list_up[[i]] <- edgeR_list_GO_up[[i]][1]
  DE_list_up[[i]] <- unlist(DE_list_up[i], use.names = FALSE)
  
  DE_list_down[[i]] <- edgeR_list_GO_down[[i]][1]
  DE_list_down[[i]] <- unlist(DE_list_down[i], use.names = FALSE)
}
#This loop runs gseGO(), enrichGO(), and groupGO() for all contrasts. Currently 
# only using edgeR outputs

for (i in 19:ncol(contrasts)){
  print(paste("Start of loop", i))
  # Rank genes for GSEA and visualization
  ranked_genes <- mutate(edgeR_list[[i]], rankstat = -log10(PValue) * sign(logFC)) |> 
    sort_by(~rankstat, decreasing = TRUE)
  geneList <- pull(ranked_genes, rankstat)
  names(geneList) <- pull(ranked_genes, genes)
  
  #GSEA
  gseGO_ALL <- gseGO(geneList = geneList,
                     OrgDb = tgu,
                     keyType = "SYMBOL",
                     ont = "ALL",
                     pvalueCutoff = 0.05, # Default 
                     verbose = FALSE
  )
  if (nrow(gseGO_ALL) > 0){
    png(filename = paste("./GSEA plots/",
                         colnames(contrasts)[i], ",All GO categories.dotplot.png",
                         sep = ""),
        height = 1100, width = 800)
    dotplot(gseGO_ALL, 
            title = paste(colnames(contrasts)[i], ", All GO categories", sep = ""), 
            split = "ONTOLOGY",
            showCategory = 15) |>
      print()
    dev.off()
    p <- cnetplot(gseGO_ALL,
                  color.params = list(foldChange = geneList),
                  circular = TRUE) + ggtitle (colnames(contrasts)[i])
    print(p)
  } else { print(paste("No gseGO() results for ", colnames(contrasts)[i]))}
  #Table of genes w/GO
  geneTable_up <- AnnotationDbi::select(tgu,
                                        keys = DE_list_up[[i]],
                                        keytype = "SYMBOL",
                                        columns = c("GENENAME","GID", "ENTREZID")
  )
  geneTable_down <- AnnotationDbi::select(tgu,
                                          keys = DE_list_down[[i]],
                                          keytype = "SYMBOL",
                                          columns = c("GENENAME","GID", "ENTREZID")
  )
  if (isEmpty(geneTable_up) && isEmpty(geneTable_down)){
    print("No DE genes founds...")
    next}
  if (!isEmpty(geneTable_up)){
    #GO over-representation analysis up
    eGO_up <- enrichGO(
      gene = geneTable_up[,1],
      OrgDb = tgu,
      keyType = "SYMBOL",
      ont = "ALL",
      pAdjustMethod = "BH",
      universe = gene_background,
      readable = TRUE,
    )
    if (!is.null(eGO_up) && nrow(eGO_up) > 0){
      png(filename = paste("./enrichGO plots/",
                           colnames(contrasts)[i], ",All GO categories up.dotplot.png",
                           sep = ""),
          height = 1100, width = 800)
      dotplot(eGO_up,
              color = "qvalue",
              title = paste(colnames(contrasts)[i], ", up", sep = "") |> print()
      )
      dev.off()
    } else print(paste("No enrichGO() up results for ", colnames(contrasts)[i]))
    #GO classification up
    gGO <- groupGO(
      gene = geneTable_up[,4][!is.na(geneTable_up[,4])],
      OrgDb = tgu,
      keyType = "ENTREZID",
      ont = "MF", #Must select one of MF, BP, or CC
      level = 2,
      readable= TRUE
    )
    if (sum(gGO[,3]) > 1) {
      p <- cnetplot(gGO,
                    #color.params = list(foldChange = geneList),
      ) + ggtitle(paste(colnames(contrasts)[i], "Molecular function, up"))
      print(p)
    } else {print("No MF up results")}
    gGO <- groupGO(
      gene = geneTable_up[,4][!is.na(geneTable_up[,4])],
      OrgDb = tgu,
      keyType = "ENTREZID",
      ont = "BP", #Must select one of MF, BP, or CC
      level = 2,
      readable= TRUE
    )
    if (sum(gGO[,3]) > 1) {
      p <- cnetplot(gGO,
                    #color.params = list(foldChange = geneList),
      ) + ggtitle(paste(colnames(contrasts)[i], "Biological process, up"))
      print(p)
    } else { print("No BP up results")}
    gGO <- groupGO(
      gene = geneTable_up[,4][!is.na(geneTable_up[,4])],
      OrgDb = tgu,
      keyType = "ENTREZID",
      ont = "CC", #Must select one of MF, BP, or CC
      level = 3,
      readable= TRUE
    )
    if (sum(gGO[,3]) > 1) {
      p <- cnetplot(gGO,
                    #color.params = list(foldChange = geneList),
      ) + ggtitle(paste(colnames(contrasts)[i], "Cellular component, up"))
      print(p)
    } else { print ("No CC up results")}
  } else { print(paste("No upregulated genes found for ", colnames(contrasts)[i]))}
  if (!isEmpty(geneTable_down)){
    #GO over-representation analysis down
    eGO_down <- enrichGO(
      gene = geneTable_down[,1],
      OrgDb = tgu,
      keyType = "SYMBOL",
      ont = "ALL",
      pAdjustMethod = "BH",
      universe = gene_background,
      readable = TRUE
    )
    if (!is.null(eGO_down)&& nrow(eGO_down) > 0){
      png(filename = paste("./enrichGO plots/",
                           colnames(contrasts)[i], ",All GO categories down.dotplot.png",
                           sep = ""),
          height = 1100, width = 800)
      dotplot(eGO_down,
              color = "qvalue",
              title = paste(colnames(contrasts)[i], ", down", sep = "") |> print()
      )
      dev.off()
    } else { print(paste("No enrichGO() results found for ", colnames(contrasts)[i], ", down"))}
    #GO classification down
    gGO <- groupGO(
      gene = geneTable_down[,4][!is.na(geneTable_down[,4])],
      OrgDb = tgu,
      keyType = "ENTREZID",
      ont = "MF", #Must select one of MF, BP, or CC
      level = 2,
      readable= TRUE
    ) 
    if (sum(gGO[,3]) > 1) {
      p <- cnetplot(gGO,
                    #color.params = list(foldChange = geneList),
      ) + ggtitle(paste(colnames(contrasts)[i], "Molecular function, down"))
      plot(p)
    } else { print("No MF down results")}
    gGO <- groupGO(
      gene = geneTable_down[,4][!is.na(geneTable_down[,4])],
      OrgDb = tgu,
      keyType = "ENTREZID",
      ont = "BP", #Must select one of MF, BP, or CC
      level = 2,
      readable= TRUE
    )
    if (sum(gGO[,3]) > 1) {
      p <- cnetplot(gGO,
                    #color.params = list(foldChange = geneList),
      ) + ggtitle(paste(colnames(contrasts)[i], "Biological process, down"))
      plot(p)
    } else { print("No BP down results")}
    gGO <- groupGO(
      gene = geneTable_down[,4][!is.na(geneTable_down[,4])],
      OrgDb = tgu,
      keyType = "ENTREZID",
      ont = "CC", #Must select one of MF, BP, or CC
      level = 3, #Higher than 3 is not very informative
      readable= TRUE
    )
    if (sum(gGO[,3]) > 1) {
      p <- cnetplot(gGO,
                    #color.params = list(foldChange = geneList),
      ) + ggtitle(paste(colnames(contrasts)[i], "Cellular component, down"))
      plot(p)
    } else {print("No CC down results")}
  } else { print(paste("No downregulated genes found for ", colnames(contrasts)[i]))}
}
