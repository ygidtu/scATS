library(e1071)

# calculate distance matrix
cal_distmat = function(data, keep_rate_cutoff=0.8){
  n_cell = dim(data$X)[1]
  keep_rate = apply(data$X>0,2,sum)/n_cell
  gene_inds = which(keep_rate>keep_rate_cutoff)
  
  # calculate distance matrix
  dist_mat = hamming.distance(data$Z[,gene_inds])
  
  return(dist_mat)
}

get_dropout_rate = function(data){
  n_cell = dim(data$X)[1]
  dropout_rate = apply(data$X==0,2,sum)/n_cell
  return(dropout_rate)
}

# process the output of odgmm with uniform component
# extract the cell labels (component id) of this gene
proc_res = function(res){
  # Z = res$Z        # n_cell x K+2 matrix,
  #                  # Note that the first is dropout and the last component in Z is uniform component
  mus = res$mus    # K+1 components, the first is dropout
  sds = res$sigmas # K+1 components, the first is dropout
  ws = res$ws      # K+2 components, the first is dropout, the last is uniform
  
  # y = apply(Z,1,which.max)
  # y = y-1      
  y = res$label - 1   # 0 means dropout, K+1 is uniform component
  
  return(list(label=y, mus=mus[-1], sds=sds[-1], ws=ws[-1] ))
}

# random impute "n_impute" missing values for gene "gene_id" using GMM component labels
rand_impute = function(data, gene_id, distmat, n_impute){
  y_impute_val = data$X[,gene_id]
  label_vec = data$Z[,gene_id]  # 0 to K+1, 0 is dropout, K+1 is uniform 
  y_mus = data$mus[[gene_id]]   # K-dimentional vector, dropout component removed
  max_x = max(data$X[,gene_id])  
  k = data$k  # k nearest neighbours
  
  dropout_inds = which(label_vec==0)
  if(n_impute>length(dropout_inds)){
    warning("Number of imputing values exceed total missing values.")
    n_impute = length(dropout_inds)
    samp_inds = dropout_inds
  }else{
    samp_inds = sort(sample(dropout_inds, n_impute))
  }
  
  valid_inds = setdiff(seq(length(label_vec)), dropout_inds)
  if(length(valid_inds)==0) return(y_impute_val)
  
  NN = distmat[samp_inds,valid_inds,drop=FALSE]
  NN = t(apply(NN, 1, order))
  
  # predicted label for the sampled indices
  if(length(valid_inds)<k){
    k=length(valid_inds)
  }
  pred = apply(NN[, 1:k, drop=FALSE], 1, function(nn){
    tmp = rle(sort(label_vec[valid_inds[nn]]))
    return(tmp$values[which.max(tmp$lengths)])  # rle(sort(c(3,4,4,2,5,7,3,3)))
  })
  
  label_vec[samp_inds] = pred
  
  K = length(y_mus)
  unif_inds = which(pred==(K+1))
  
  if(length(unif_inds)>0){
    oth_inds = which(pred<=K)
    y_impute_val[samp_inds[unif_inds]] = runif(length(unif_inds), 0, max_x)
    y_impute_val[samp_inds[oth_inds]] =  y_mus[pred[oth_inds]]
  }else{
    y_impute_val[samp_inds] = y_mus[pred]
  }
  
  return(y_impute_val)
}

# random impute "n_impute" missing values for gene "gene_id" using GMM component labels
rand_impute_v1 = function(data, gene_id, distmat, n_impute){
  y_impute_val = data$X[,gene_id]
  label_vec = data$Z[,gene_id]  # 0 to K+1, 0 is dropout, K+1 is uniform 
  y_mus = data$mus[[gene_id]]   # K-dimentional vector, dropout component removed
  max_x = max(data$X[,gene_id])  
  k = data$k  # k nearest neighbours
  
  dropout_inds = which(label_vec==0)
  if(n_impute>length(dropout_inds)){
    warning("Number of imputing values exceed total missing values.")
    n_impute = length(dropout_inds)
    samp_inds = dropout_inds
  }else{
    samp_inds = sort(sample(dropout_inds, n_impute))
  }
  
  NN = distmat[samp_inds, ,drop=FALSE]
  NN = t(apply(NN, 1, order))
  NN = NN[,1:k,drop=FALSE]
  
  pred = apply(NN, 1, function(nn){
    tmp = rle(sort(label_vec[nn]))
    return(tmp$values[which.max(tmp$lengths)])  # rle(sort(c(3,4,4,2,5,7,3,3)))
  })
  
  label_vec[samp_inds] = pred
  
  valid_inds = which(pred!=0)
  if(length(valid_inds)==0){
    return(y_impute_val)
  }else{
    samp_inds = samp_inds[valid_inds]
    pred = pred[valid_inds]
  }
  
  K = length(y_mus)
  unif_inds = which(pred==(K+1))
  
  if(length(unif_inds)>0){
    oth_inds = which(pred<=K)
    y_impute_val[samp_inds[unif_inds]] = runif(length(unif_inds), 0, max_x)
    y_impute_val[samp_inds[oth_inds]] =  y_mus[pred[oth_inds]]
  }else{
    y_impute_val[samp_inds] = y_mus[pred]
  }
  
  return(y_impute_val)
}

# calculate the pdf of the gmm+uni mixture model at x
# unif_max: uniform compoent
gmm_uni_pdf = function(x, mus, sigmas, ws, unif_max, log=FALSE){
  K = length(mus)
  stopifnot(length(ws)==(K+1))   # last component is uniform component
  N = length(x)
  y = rep(0,N)
  for(k in seq(K)){
    y = y + ws[k] * dnorm(x,mus[k],sigmas[k])
  }
  y = y + ws[K+1]/unif_max*(x<unif_max)
  
  if(log){
    return(log(y))
  }else{
    return (y)
  }
}

# extract useful information from data for gene_id
ex_para = function(data, gene_id){
  ws = data$ws[[gene_id]][-1]  # first component is dropout
  mus = data$mus[[gene_id]][-1]
  sds = data$sds[[gene_id]][-1]
  max_x = max(data$X[,gene_id])
  return(list(ws=ws/sum(ws), mus=mus, sds=sds, max_x=max_x))
}

# calculate the expectation of the gmm+uniform mixture
cal_gmm_uni_exp = function(mus, ws, unif_expectation){
  K = length(mus)
  return( sum(mus*ws[1:K] + unif_expectation*ws[K+1]) )
}

# calculate the KL divergence between para_l and para_r, by shifting para1 "delta_x" distance to the right 
# para_l: parameters for the gmm that is on the left, i.e. lower expectation
# para_r: parameters for the gmm on the right, i.e. high expectation
# x: x for calculating the pdf
cal_shift_KL = function(para_l, para_r, x, delta_x_arr){
  yr = gmm_uni_pdf(x, para_r$mus, para_r$sds, para_r$ws, para_r$max_x)
  yr = yr + 1e-6
  
  len = length(delta_x_arr)
  kl_arr = rep(0,len)
  
  for(i in seq(len)){
    delta_x = delta_x_arr[i]
    yl = gmm_uni_pdf(x, para_l$mus+delta_x, para_l$sds, para_l$ws, para_l$max_x+delta_x)
    yl = yl + 1e-6
    kl_arr[i] = sum((log(yl)-log(yr))*yl)
  }
  return( delta_x_arr[which.max(kl_arr)] )
}

# calculate the profile shift for gene_id
# positive value means shift data1, negative means shift data2
cal_profile_shift = function(data1, data2, gene_id){
  para1 = ex_para(data1,gene_id)
  para2 = ex_para(data2,gene_id)
  exp1 = cal_gmm_uni_exp(para1$mus, para1$ws, para1$max_x/2)
  exp2 = cal_gmm_uni_exp(para2$mus, para2$ws, para2$max_x/2)
  
  max_x = max(para1$max_x, para2$max_x)
  x = seq(0,max_x,1000)
  exp_diff = abs(exp1-exp2)*2
  
  delta_x_arr = seq(0,exp_diff,0.1)
  if(length(delta_x_arr)==1 && delta_x_arr==0){
    return(0)
  }
    
  if(exp1<exp2){
    return(cal_shift_KL(para1, para2, x, delta_x_arr))
  }else{
    return(-cal_shift_KL(para2, para1, x, delta_x_arr))
  }
}

gen_data = function(GMM_res, gene_expr_mat){
  stopifnot( dim(GMM_res[[1]]$Z)[2] == (length(GMM_res[[1]]$mus)+1)  )  #  GMM with uniform component
  
  # GMM_res is a list, 1 x n_gene
  n_gene = length(GMM_res)
  n_cell = dim(gene_expr_mat)[1]
  stopifnot(dim(gene_expr_mat)[2]==n_gene)
  stopifnot(dim(GMM_res[[1]]$Z)[2] == (length(GMM_res[[1]]$mus)+1)  )  #  GMM with uniform component
  
  Z = matrix(0, nrow = n_cell, ncol = n_gene)   # gmm label matrix
  mus = vector("list", length = n_gene)
  sds = vector("list", length = n_gene)
  ws = vector("list", length = n_gene)
  
  ##GMM_res is a list contain res 
  for(i in seq(n_gene)){
    cat("proc gene",i,"\n")
    tmpres = proc_res(GMM_res[[i]])
    Z[,i] = tmpres$label
    mus[[i]] = tmpres$mus  # K-dimensional component
    sds[[i]] = tmpres$sds
    ws[[i]]  = tmpres$ws
  }
  
  data = list(n_cell = n_cell,
              n_gene = n_gene,
              X = gene_expr_mat,   # n_cell x n_gene matrix, log(x+1) where x is the raw gene count of the cell
              Z = Z,   # n_cell x n_gene matrix, each element is the cell's label (GMM component) of the gene
              mus=mus, # 1 x n_gene list, each contains the GMM means of each gene, K-dimentional vector (no dropout component)
              sds=sds, # 1 x n_gene list, each contains the GMM sds of each gene, K-dimentional vector
              ws =ws, # 1 x n_gene list, each contains the GMM weights of each gene, K+1-dimentional vector
              k=10         # k nearest label
  )
  
  return(data)
}

###############  main code ##################################

# some irrelevant code to fix bug
fix_res = function(res_list, gene_expr_mat){   # n_cell x n_gene
  flag_vec = as.logical( lapply(res_list, FUN = function(x){!is.null(x$Z)}) )
  for(i in which(flag_vec)){
    tmp = res_list[[i]]
    z_ind = which(names(tmp) == "Z")
    tmp = tmp[-z_ind]
    tmp$label = rep(1,length(gene_expr_mat[,i]))
    res_list[[i]] = tmp
  }
  return(res_list)
}

###############  step 1: prepare data ##################################
# GMM_res1 = GMM_res_batch1_keeprate_0.1  # a list, 1 x n_gene
# gene_expr_mat1 = unname( as.matrix(t(X_tab_batch1_log_keeprate_0.1)) ) # n_cell x n_gene
# GMM_res1 = fix_res(GMM_res1, gene_expr_mat1)
# 
# GMM_res2 = GMM_res_batch2_keeprate_0.1
# gene_expr_mat2 = unname( as.matrix(t(X_tab_batch2_log_keeprate_0.1)) )
# GMM_res2 = fix_res(GMM_res2, gene_expr_mat2)

# data1 = gen_data(GMM_res1, gene_expr_mat1)
# data2 = gen_data(GMM_res2, gene_expr_mat2)

stopifnot(data1$n_gene==data2$n_gene)

# calculate distance matrix
cat("start calculating distance matrix\n")
keep_rate_cutoff = 0.3    # gene index used to calculate distance matrix
# distmat1 = cal_distmat(data1, keep_rate_cutoff)
# distmat2 = cal_distmat(data2, keep_rate_cutoff)
cat("distance matrix ready\n")

###############  step 2: impute dropout  ##################################
dropout_rate_vec1 = get_dropout_rate(data1)   # all genes
dropout_rate_vec2 = get_dropout_rate(data2)
keep_rate_vec1 = 1 - dropout_rate_vec1
keep_rate_vec2 = 1 - dropout_rate_vec2
tmpdiff = dropout_rate_vec1 - dropout_rate_vec2
tmp1 = round(abs(data1$n_cell*tmpdiff))
tmp2 = round(abs(data2$n_cell*tmpdiff))
n_impute_arr = rep(0,length(dropout_rate_vec1))
n_impute_arr[tmpdiff>0] = tmp1[tmpdiff>0]
n_impute_arr[tmpdiff<0] = tmp2[tmpdiff<0]

Z1 = data1$Z
Z2 = data2$Z
X1 = data1$X
X2 = data2$X
for(i in seq(length(dropout_rate_vec1))){
  if(i%%100==0){
    cat(paste0("impute ",i,"th gene.\n"))
  }
  if(n_impute_arr[i]<2) next()
  # if(dropout_rate_vec1[i]>0.9 || dropout_rate_vec2[i]>0.9) next()
  if(dropout_rate_vec1[i]>dropout_rate_vec2[i]){
    # impute gene i for batch 1
    # X1[,i] = rand_impute(data1, i, distmat1, n_impute_arr[i])
    X1[,i] = rand_impute_v1(data1, i, distmat1, n_impute_arr[i])
  }else{
    # impute gene i for batch 2
    # X2[,i] = rand_impute(data2, i, distmat2, n_impute_arr[i])
    X2[,i] = rand_impute_v1(data2, i, distmat2, n_impute_arr[i])
  }
}

###############  step 3: profile matching  ##################################
stopifnot(data1$n_gene==data2$n_gene)
n_gene = data1$n_gene
# valid_gene_inds = which(keep_rate_vec1>0.2 & keep_rate_vec2>0.2)
# valid_gene_inds = seq(length(dropout_rate_vec1))
for(i in seq(length(dropout_rate_vec1))){
  tmp_shift = cal_profile_shift(data1, data2, i)
  if(tmp_shift>=0){
    X1[,i] = X1[,i] + tmp_shift
  }else{
    X2[,i] = X2[,i] + abs(tmp_shift)
  }
}
# X1 = X1[,valid_gene_inds]
# X2 = X2[,valid_gene_inds]

###############  step 4: plot PCA  ##################################
# mtcars.pca <- prcomp(mtcars[,c(1:7,10,11)], center = TRUE,scale. = TRUE)
# summary(mtcars.pca)
library(ggfortify)

valid_gene_inds = which(abs(keep_rate_vec1-keep_rate_vec2)<0.1)

X = rbind(X1[,valid_gene_inds],X2[,valid_gene_inds])
lab_vec = matrix(c(rep('1',dim(X1)[1]), rep('2',dim(X2)[1])),ncol=1)
# colnames(lab_vec) = "lab"
# autoplot(prcomp(X),lab_vec, colour="lab")

library(ggplot2)
library(Rtsne)
tsne_out <- Rtsne(X) # Run TSNE
d_tsne = as.data.frame(tsne_out$Y)
ggplot(d_tsne, aes(x=V1, y=V2, color=lab_vec)) +
  geom_point(size=0.1)

cov_mat_inv = solve(cov(X))
XO = X%*%cov_mat_inv
tsne_out1 <- Rtsne(XO) # Run TSNE
d_tsne = as.data.frame(tsne_out1$Y)
ggplot(d_tsne, aes(x=V1, y=V2, color=lab_vec)) +
  geom_point(size=0.1)

