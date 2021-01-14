dexde_test <- function(obj,
                       idents.1,
                       idents.2 = NULL,
                       annot,
                       cores=1,
                       assay = 'RNA',
                       num.splits = 6) {
  set.seed(1)
  
  pa_use <- annot$pa_site
  
  cells.1 = names(Seurat::Idents(obj))[which(Seurat::Idents(obj) == idents.1)]
  if (is.null(idents.2)) {
    cells.2 = names(Seurat::Idents(obj))[which(Seurat::Idents(obj) != idents.1)]
  } else {
    cells.2 = names(Seurat::Idents(obj))[which(Seurat::Idents(obj) == idents.2)]
  }
  
  cells.sub.1 = split(cells.1, sort(1:length(cells.1) %% num.splits))
  data_to_test.1 = matrix(, nrow = length(pa_use), ncol = length(cells.sub.1))
  
  cells.sub.2 = split(cells.2, sort(1:length(cells.2) %% num.splits))
  data_to_test.2 = matrix(, nrow = length(pa_use), ncol = length(cells.sub.2))
  
  
  for (i in 1:length(cells.sub.1)) {
    this.set <- cells.sub.1[[i]]
    sub.matrix <-
      Seurat::GetAssayData(obj, slot = "counts", assay = assay)[pa_use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x)
        sum(x)))
      data_to_test.1[, i] <- this.profile
    } else {
      data_to_test.1[, i] <- sub.matrix
    }
  }
  
  rownames(data_to_test.1) <- pa_use
  colnames(data_to_test.1) <-
    paste0("Population1_", 1:length(cells.sub.1))
  
  
  
  
  
  
  for (i in 1:length(cells.sub.2)) {
    this.set <- cells.sub.2[[i]]
    sub.matrix <-
      Seurat::GetAssayData(obj, slot = "counts", assay = assay)[pa_use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x)
        sum(x)))
      data_to_test.2[, i] <- this.profile
    } else {
      data_to_test.2[, i] <- sub.matrix
    }
  }
  
  rownames(data_to_test.2) <- pa_use
  colnames(data_to_test.2) <-
    paste0("Population2_", 1:length(cells.sub.1))
  
  
  peak.matrix <- cbind(data_to_test.1, data_to_test.2)
  rownames(peak.matrix) <-
    stringr::str_replace_all(rownames(peak.matrix), ':', '_')
  
  sampleTable <-
    data.frame(row.names = c(colnames(data_to_test.1), colnames(data_to_test.2)),
               condition = c(rep("target", ncol(data_to_test.1)),
                             rep("comparison", ncol(data_to_test.2))))
  
  exon_info <- annot[, c('pa_site', 'annot.symbol')]
  exon_info$pa_site <-
    stringr::str_replace_all(exon_info$pa_site, ':', '_')
  
  library(DEXSeq)
  library(magrittr)
  
  dxd <- DEXSeqDataSet(
    peak.matrix[, rownames(sampleTable)],
    sampleTable,
    design = ~ sample + exon + condition:exon,
    exon_info$pa_site,
    exon_info$annot.symbol,
    featureRanges = NULL,
    transcripts = NULL,
    alternativeCountData = NULL
  )
  
  require(magrittr)
  if (cores != 1) {
    BPPARAM = BiocParallel::MulticoreParam(cores)
    dxd %<>%
      estimateSizeFactors %<>%
      estimateDispersions(BPPARAM = BPPARAM) %<>%
      testForDEU(BPPARAM = BPPARAM) %<>%
      estimateExonFoldChanges(BPPARAM = BPPARAM)

  } else {
    dxd %<>%
      estimateSizeFactors %<>%
      estimateDispersions %<>%
      testForDEU %<>%
      estimateExonFoldChanges
  }
  
  dxr1 = DEXSeq::DEXSeqResults(dxd)
  dxr1 <-
    dxr1[, c("groupID",
             "exonBaseMean",
             "pvalue",
             "padj",
             "log2fold_target_comparison")]
  dxr1 <- as.data.frame(dxr1)
  pa_use <- stringr::str_replace_all(sapply(strsplit(rownames(dxr1), split = ":"), "[[", 2),
                                     "_", ":")
  dxr1[["pct.1"]] <-
    apply(
      X = GetAssayData(obj, "counts", assay = assay)[pa_use, Idents(obj) == idents.1],
      MARGIN = 1,
      FUN = Seurat:::PercentAbove,
      threshold = 0
    )
  if (is.null(idents.2)) {
    dxr1[["pct.2"]] <-
      apply(
        X = GetAssayData(obj, "counts", assay = assay)[pa_use, Idents(obj) != idents.1],
        MARGIN = 1,
        FUN = Seurat:::PercentAbove,
        threshold = 0
      )
    dxr1[["cluster"]] <- glue::glue('{idents.1}_vs_Other')
  } else {
    dxr1[["pct.2"]] <-
      apply(
        X = GetAssayData(obj, "counts", assay = assay)[pa_use, Idents(obj) == idents.2],
        MARGIN = 1,
        FUN = Seurat:::PercentAbove,
        threshold = 0
      )
    dxr1[["cluster"]] <- glue::glue('{idents.1}_vs_{idents.2}')
  }
  
    
  return(dxr1)
}
