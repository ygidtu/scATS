arg <- commandArgs(T)
file <- arg[1]
source('scAPA.R')

library(doMC)

registerDoMC(10)

# file <- "../Genes_new/10_77268681_77269575.txt"
file = '/mnt/raid61/Personal_data/zhangyiming/code/afe/tests/extract/test.txt'
set.seed(1)


df <- read.table(file, row.names = NULL,  header = T)
# df$V7 <- df$V1 + df$V2 -1

# df$V7[df$V5 == "no"] <- NA
# df$V7 <- df$V7 + abs(min(df$V1)) + 1
# df$V1 <- df$V1 + abs(min(df$V1)) + 1


# L <- df[which.max(df$V1), 'V1'] + df[which.max(df$V2), 'V2'] + 50

# L2 <- max(df$V7[!is.na(df$V7)])
# L <- max(L1, L2)

# df = df[1:10000, ]

# df = df[df$ReadStart > df$UTRStart & df$ReadStart < df$UTREnd, ]

df <- df[df$UTRStart == "2841486" & df$UTREnd == "2841973",]
df <- df[df$StartLocInUTR > -500,]
df$LenPA <- NA

r1_utr_st_arr=df$StartLocInUTR # start location of each read on UTR part

r1_utr_st_arr <- r1_utr_st_arr
r1_len_arr=df$LenInUTR # length of each read on UTR part
polya_len_arr=df$LenPA # length of polyA
# pa_site_arr=df$PASite, # pa site locations of reads
L = mean(df$UTREnd - df$UTRStart)

r1_utr_st_arr <- r1_utr_st_arr + L

utr_l <- max(r1_utr_st_arr) + max(r1_len_arr) + 300

source('/mnt/raid61/Personal_data/zhangyiming/code/afe/jl/scAPA.R')
res = apamix(
  n_max_apa = 5,
  n_min_apa = 1,
  r1_utr_st_arr,
  r1_len_arr,
  polya_len_arr,
  polya_len_arr,
  polya_len_arr,
  L = utr_l,
  mu_f=270,
  sigma_f = 30,
  min_ws = 0.01,
  debug=T
)


res = apamix(
  n_max_apa = 5,
  n_min_apa = 1,
  r1_utr_st_arr=df$V1, # start location of each read on UTR part
  r1_len_arr=df$V2, # length of each read on UTR part
  r2_len_arr=df$V3, # length of each read on polyA part
  polya_len_arr=df$V4, # length of polyA
  # pa_site_arr=df$PASite, # pa site locations of reads
  L = df[which.max(df$V1), 'V1'] + df[which.max(df$V2), 'V2'] + 50,
  mu_f=270,
  sigma_f = 30,
  min_ws = 0.01,
  debug=T
)


label <- gsub(".txt","",basename(file))
saveRDS(res, file = glue::glue('extract/{label}.Rds'))
# res <- readRDS(glue::glue('{label}.Rds'))

x = df$V1
l = df$V2
r = df$V3
s = df$V6
L = L
LA = 120
n_frag = dim(df)[1]

s[is.na(s)]=60
coverage_cnt = rep(0,L+LA)
png(glue::glue('extract/{label}.png'), res=600, height = 6, width = 12, units = "in")
for(i in seq(n_frag)){
  coverage_cnt[x[i]:(x[i]+l[i]-1)] = coverage_cnt[x[i]:(x[i]+l[i]-1)] + 1
  coverage_cnt[(L+s[i]-r[i]+1):(L+s[i])] = coverage_cnt[(L+s[i]-r[i]+1):(L+s[i])] + 1
}

plot(seq(L+LA),coverage_cnt,type="s",
     # xlab = paste("UTR | polyA","ws=", paste(sprintf("%.2f",ws), collapse=" ")), 
     xlab = paste0("UTR | polyA"), 
     ylab = "coverage count",
     main = paste0(label, " ", paste( round(res$ws,digits = 2), sep = " ", collapse = " ") ) ,
     xaxt='n')
abline(v=L, col="green",lwd=3)

ymax = par("usr")[4]

for(i in seq(length(res$alpha_arr)) ){
  st = res$alpha_arr[i] - res$beta_arr[i]
  en = res$alpha_arr[i] + res$beta_arr[i]
  x = c(st,st,en,en)
  y = c(0,ymax,ymax,0)
  lines(x,y,type='S',col='blue',lwd=2,lty=2)
  abline(v=res$alpha_arr[i], col="yellow",lwd=2,lty=1)
}

axis(side = 1, at = c(seq(0,L+LA,100)), las=2 )

jrt = as.data.frame(table(df$V7[!is.na(df$V7)]))
loc = as.numeric(levels(jrt[,1]))
# pmf = jrt[,2]/max(jrt[,2])
lines(loc, jrt[,2], type='h', col='red', lwd=1)

dev.off()

