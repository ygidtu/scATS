
library(Seurat)
library(dplyr)

dat <- readRDS('/mnt/raid64/ATS/Personal/zhangyiming/cage_sel/200/seurat_obj.rds')

annot <- as.data.frame(readRDS('annot.Rds'))

annot[['pa_site']] <- gsub('chr','',paste(annot$seqnames, paste(annot$start, annot$end, sep = '-'), annot$strand, sep = ':'))
rownames(annot) <- annot[['pa_site']]
annot[['annot.symbol']] <- annot[['geneID']]
annot <- annot[rownames(dat),]

annot %>%
  subset(!is.na(geneID)) %>%
  group_by(geneID) %>%
  filter(n() > 1) %>% as.data.frame -> annot_sub

source('/mnt/raid61/Personal_data/zhangyiming/code/afe/cage/de_utils.R')
Idents(dat) <- 'CellType'


cell_type_res = list()
for (i in unique(dat@meta.data$CellType)) {
  print(i)

  cell_type_res[[i]] = dexde_test(dat,
                       idents.1 = i,
                       # idents.2 = 'erythroid',
                       idents.2 = NULL,
                       annot_sub,
                       cores=1,
                       assay = 'RNA',
                       num.splits = 6)
}

saveRDS(cell_type_res, "cell_type_dex.rds")

dat@meta.data$temp_group = paste(dat@meta.data$CellType, dat@meta.data$group, sep = "|")
Idents(dat) <- 'temp_group'
res = list()
for (i in unique(dat@meta.data$CellType)) {
    print(i)

    for (j in c("NCovM", "NCovS")) {
        tryCatch({
          res[[paste(i, paste(j, "HN", sep = " - "), sep  = "|")]] = dexde_test(
            dat,
            idents.1 = paste(i, j, sep = "|"),
            idents.2 = paste(i, "HN", sep = "|"),
            annot_sub,
            cores=1,
            assay = 'RNA',
            num.splits = 6
          )
        }, error = function(e) {})
    }
    
    for (j in c("PCovM", "PCovS")) {
        tryCatch({
          res[[paste(i, paste(j, "HP", sep = " - "), sep  = "|")]] = dexde_test(
            dat,
            idents.1 = paste(i, j, sep = "|"),
            idents.2 = paste(i, "HP", sep = "|"),
            annot_sub,
            cores=1,
            assay = 'RNA',
            num.splits = 6
          )
        }, error = function(e) {})
    }

    tryCatch({
      res[[paste(i, "NCovS - NCovM", sep  = "|")]] = dexde_test(
        dat,
        idents.1 = paste(i, "NCovS", sep = "|"), 
        idents.2 = paste(i, "NCovM", sep = "|"),
        annot_sub,
        cores=1,
        assay = 'RNA',
        num.splits = 6
      )
    }, error = function(e) {})

    tryCatch({
      res[[paste(i, "PCovS - PCovM", sep  = "|")]] = dexde_test(
        dat,
        idents.1 = paste(i, "PCovS", sep = "|"), 
        idents.2 = paste(i, "PCovM", sep = "|"),
        annot_sub,
        cores=1,
        assay = 'RNA',
        num.splits = 6
      )
    }, error = function(e) {})
}

saveRDS(res, "groups_dex.rds")