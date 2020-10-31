# simulate data for APA

############  auxiliary function  ######################
# generate integers that are normally distributed
gen_frag_len = function(n, mu=0, sigma=1){
  mu = round(mu)
  sigma = round(sigma)
  ymin = mu - 3*sigma
  ymax = mu + 3*sigma
  
  y = round(rnorm(n, mu, sigma))
  y[y<ymin] = ymin
  y[y>ymax] = ymax
  
  return(y)
}

############  part 1:  parameter configurations  ######################
utr_len = 1000
polya_len = 250

min_polya_len = 50
min_pa_site = 150

n_frag = 1000

mu_frag = 300
sig_frag = 20

mu_r1 = 100
sig_r1 = 10

mu_r2 = 30
sig_r2 = 5

stopifnot(min_polya_len>( mu_r2 + 3*sig_r2))
stopifnot(min_pa_site>( mu_r1 + 3*sig_r1 ))


############  part 2:  PA sites  ######################
n_pa_site = 2
min_gap = 100*rep(1,n_pa_site)  # minimum distance between PA sites
min_gap[1] = 150

gap_size = sum(min_gap)
stopifnot(gap_size<utr_len)

theta_arr = sort(sample(utr_len-gap_size, n_pa_site))
theta_arr = theta_arr + cumsum(min_gap)

ws = 1+runif(n_pa_site)
ws = ws/sum(ws)

############  part 3:  fragments and reads ######################

reads_boarder = round(n_frag * ws)
reads_boarder = cumsum(reads_boarder)
reads_boarder[length(reads_boarder)] = n_frag
reads_boarder = c(0,reads_boarder)

frag_label = rep(0,n_frag)

polya_part_len_arr = min_polya_len + sample(polya_len-min_polya_len,n_frag,replace = T)   # s
frag_len_arr = gen_frag_len(n_frag, mu_frag, sig_frag)  # theta - x + 1 + s
read_utr_st_arr = rep(0,n_frag)    # x

theta_vec = rep(0,n_frag)
for(i in seq(n_pa_site)){
  tmpinds = seq((reads_boarder[i]+1),reads_boarder[i+1])
  theta_vec[tmpinds] = theta_arr[i]
  frag_label[tmpinds] = i
  read_utr_st_arr[tmpinds] = -frag_len_arr[tmpinds] + theta_arr[i] + 1 + polya_part_len_arr[tmpinds]
}

r1_len_arr = gen_frag_len(n_frag, mu_r1, sig_r1)   # l
r2_len_arr = gen_frag_len(n_frag, mu_r2, sig_r2)   # r

valid_ind1 = ( r1_len_arr <= (theta_vec - read_utr_st_arr + 1) )
valid_ind2 = ( r2_len_arr <= polya_part_len_arr )
valid_ind3 = ( read_utr_st_arr > 0 )

valid_inds = (valid_ind1 & valid_ind2 & valid_ind3)

n_frag = sum(valid_inds)
theta = theta_arr
z = frag_label[valid_inds]
ws = as.data.frame(table(z))$Freq
ws = ws/sum(ws)

L = utr_len
LA = polya_len
x = read_utr_st_arr[valid_inds]
theta_vec = theta_vec
l = r1_len_arr[valid_inds]
r = r2_len_arr[valid_inds]
s = polya_part_len_arr[valid_inds]

n_pre = 2
pre_inds = sort(sample(n_pa_site,n_pre))
pre_theta_arr = theta_arr[pre_inds]
pre_theta_cnt_arr = rep(0,n_pre)
tmp_arr = ( r1_len_arr >= (theta_vec - read_utr_st_arr + 1) )
for(i in seq(n_pre)){
  pre_theta_cnt_arr[i] = sum(tmp_arr[theta_vec==pre_theta_arr[i]])
}
cat('pre_theta_arr=',pre_theta_arr,'\n')
cat('pre_theta_cnt_arr=',pre_theta_cnt_arr,'\n')


############  part 4:  plot simulated data ######################

coverage_cnt = rep(0,L+LA)
for(i in seq(n_frag)){
  coverage_cnt[x[i]:(x[i]+l[i]-1)] = coverage_cnt[x[i]:(x[i]+l[i]-1)] + 1
  coverage_cnt[(L+s[i]-r[i]+1):(L+s[i])] = coverage_cnt[(L+s[i]-r[i]+1):(L+s[i])] + 1
}

plot(seq(L+LA),coverage_cnt,type="s",
     xlab = paste("UTR | polyA","ws=", paste(sprintf("%.2f",ws), collapse=" ")), 
     ylab = "coverage count",
     main = "red: PA sites, blue: UTR&polyA boarder")
abline(v=L, col="blue",lwd=3)
for(i in theta_arr){
  abline(v=i, col="red",lwd=2,lty=2)
}
