load <- function() {
  library(CAGEfightR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(knitr)
  library(viridis)
  library(magrittr)
  library(tidyverse)
}

# BiocManager::install("CAGEfightR")
suppressPackageStartupMessages(load())


bsg <- BSgenome.Hsapiens.UCSC.hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
odb <- org.Hs.eg.db

# Necessary function
convert_ctss_to_bigwig = function(ctss_file, genomeInfo) {
  tpm_bed = rtracklayer::import(ctss_file)
  chrom_id = levels(seqnames(tpm_bed))
  if (length(grep('chr', tolower(chrom_id))) != 0) {
    chrom_id = chrom_id[nchar(chrom_id) < 8]
  } else {
    chrom_id = chrom_id[nchar(chrom_id) < 4]
  }
  tpm_bed = tpm_bed[seqnames(tpm_bed) %in% chrom_id]
  tpm_bed = GenomicRanges::GRanges(tpm_bed , seqinfo = genomeInfo)
  tpm_bed_plus = tpm_bed[tpm_bed@strand == "+",]
  tpm_bed_minus = tpm_bed[tpm_bed@strand == "-",]
  plus_id = paste(ctss_file , ".plus.bw" , sep = "")
  rtracklayer::export(object = tpm_bed_plus ,
                      plus_id,
                      format = "BigWig")
  minus_id = paste(ctss_file , ".minus.bw" , sep = "")
  rtracklayer::export(object = tpm_bed_minus,
                      minus_id,
                      format = "BigWig")
  return(c(plus_id, minus_id))
}

args = commandArgs(trailingOnly=TRUE)
search_dir <- args[1]
# search_dir = "/mnt/raid64/ATS/Personal/zhangyiming/cage_sel"
setwd(search_dir)
genome_info = 'hg38'

ctss_files <-
  list.files(path = search_dir,
             pattern = 'R1.bed$',
             full.names = T, recursive=F)

names(ctss_files) <-
  gsub(basename(ctss_files),
       pattern = '_R1.bed',
       replacement = '')


genome_cache = glue::glue('{genome_info}.Rds')

if (file.exists(genome_cache)) {
  genomeInfo <- readRDS(genome_cache)
} else {
  genomeInfo <- SeqinfoForUCSCGenome(genome_info)
  saveRDS(genomeInfo, file = genome_cache)
}

if(file.exists("CTSSs.Rds")) {
  CTSSs = readRDS("CTSSs.Rds")
} else {
  file_lst = lapply(ctss_files, function(x) {
    convert_ctss_to_bigwig(x, genomeInfo)
  })

  sample_info =
    data.frame(
      Name = names(file_lst),
      BigWigPlus = sapply(file_lst, '[[', 1),
      BigWigMinus = sapply(file_lst, '[[', 2),
      stringsAsFactors = F
    )

  bw_plus <- BigWigFileList(sample_info$BigWigPlus)
  bw_minus <- BigWigFileList(sample_info$BigWigMinus)
  names(bw_plus) <- sample_info$Name
  names(bw_minus) <- sample_info$Name

  CTSSs <- quantifyCTSSs(
    plusStrand = bw_plus,
    minusStrand = bw_minus,
    design = sample_info,
    genome = genomeInfo
  )

  CTSSs <- CTSSs %>% calcTPM() %>% calcPooled()
  saveRDS(CTSSs, file = 'CTSSs.Rds')
}

CTSSs <- readRDS("CTSSs.Rds")

BCs <- quickEnhancers(CTSSs)
BCs <- assignTxType(BCs, txModels = txdb, swap="thick")
saveRDS(BCs, file = 'BCs.Rds')

for (i in seq(0.1, 0.6, 0.1)) {
  for (j in c(1, seq(10, 100, 10))) {
    TSSs <-
      quickTSSs(CTSSs) %>%
      calcTPM() %>%
      subsetBySupport(inputAssay = "counts",
                      unexpressed = j,
                      minSamples = round(length(ctss_files) * i))

    TSSs <- assignTxID(TSSs,
                      txModels = txdb, swap = "thick")

    peak_region <- rownames(as.data.frame(rowRanges(TSSs)))
    peak_region <-
      gsub(peak_region, 
          pattern = '[:]|[;]', 
          replacement = '\t')

    # length(peak_region)

    peak_region <- gsub(peak_region, pattern = 'chr', replacement = '')

    write.table(
      peak_region,
      file = paste('peak_region', i, j,'.txt', sep = "_"),
      quote = F,
      row.names = F,
      col.names = F
    )
  }
}
TSSs <-
  quickTSSs(CTSSs) %>%
  calcTPM() %>%
  subsetBySupport(inputAssay = "counts",
                  unexpressed = 100,
                  minSamples = round(length(ctss_files) * i))

TSSs <- assignTxID(TSSs,
                   txModels = txdb, swap = "thick")

# saveRDS(TSSs, file = 'TSSs.Rds')

peak_region <- rownames(as.data.frame(rowRanges(TSSs)))
peak_region <-
  gsub(peak_region, 
       pattern = '[:]|[;]', 
       replacement = '\t')

length(peak_region)

peak_region <- gsub(peak_region, pattern = 'chr', replacement = '')

write.table(
  peak_region,
  file = 'peak_region.txt',
  quote = F,
  row.names = F,
  col.names = F
)
