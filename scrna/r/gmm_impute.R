library(e1071)

# impute missing values
# orig_data: n_cell x n_gene matrix, original (transformed) expression value 
# gmm_label_mat: n_cell x n_gene matrix, each element is the cell's gmm component label of each gene, dropout marked as 0
# gmm_mus_list: 1 x n_gene list, each element are the means for the non-dropout gmm components and the last one is the maximal expression value for the uniform component
# k: number of nearest neighbours for imputation
# return: impute_label_mat: n_cell x n_gene matrix, with dropout labels imputed
#         impute_val_mat: n_cell x n_gene matrix, with dropout values imputed. Note original values are kept as it is 
gmm_impute = function(orig_data, gmm_label_mat, gmm_mus_list, k=5){
  n_gene = dim(orig_data)[2]
  IX = gmm_label_mat
  IX_val = orig_data
  
  distmat = hamming.distance(gmm_label_mat)  # hamming distance between the rows
  
  for(i in seq(n_gene)){
    cat('Impute column ',i,"\n")
    res = gmm_impute0(distmat, orig_data[,i], gmm_label_mat[,i], gmm_mus_list[[i]], k)
    IX[,i] = res[,1]
    IX_val[,i] = res[,2]
  }
  return(list(impute_label_mat=IX,impute_val_mat=IX_val))
}

# impute missing values, may not impute if nearest neighbours all dropout
# orig_data: n_cell x n_gene matrix, original (transformed) expression value 
# gmm_label_mat: n_cell x n_gene matrix, each element is the cell's gmm component label of each gene, dropout marked as 0
# gmm_mus_list: 1 x n_gene list, each element are the means for the non-dropout gmm components and the last one is the maximal expression value for the uniform component
# k: number of nearest neighbours for imputation
# return: impute_label_mat: n_cell x n_gene matrix, with dropout labels imputed
#         impute_val_mat: n_cell x n_gene matrix, with dropout values imputed. Note original values are kept as it is 
gmm_impute_v1 = function(orig_data, gmm_label_mat, gmm_mus_list, k=10){
  n_gene = dim(orig_data)[2]
  IX = gmm_label_mat
  IX_val = orig_data
  
  distmat = hamming.distance(gmm_label_mat)  # hamming distance between the rows
  NN = t(apply(distmat, 1, order))
  
  for(i in seq(n_gene)){
    cat('Impute column ',i,"\n")
    res = gmm_impute1(NN, orig_data[,i], gmm_label_mat[,i], gmm_mus_list[[i]], k)
    IX[,i] = res[,1]
    IX_val[,i] = res[,2]
  }
  return(list(impute_label_mat=IX,impute_val_mat=IX_val))
}

# impute all missing values indicated by orig_data_vec==0, based on similarity given by distmat
# distmat: n x n matrix, n cells, distance matrix between cells
# orig_data_vec: n x 1 vector, original (transformed) expression values for this gene
# gmm_label_vec: n x 1 vector, labels (component id) for this gene, 0 is dropout, others are 1 to K+1
# gmm_mus: K dimentional vector, each element is a component mean, dropout component (first) is removed, the last of which is the largest expression value for uniform component 
# k: k nearest neighbors
# return: y with all dropout labels (0) imputed, y_impute_val is imputed missing values (original values are copied)
gmm_impute0 = function(distmat, orig_data_vec, gmm_label_vec, gmm_mus, k){
  y = gmm_label_vec
  y_mus = gmm_mus
  
  max_gene_expr = max(orig_data_vec)
    
  test_inds = which(y==0)
  train_inds = which(y!=0)
  #y_impute_val = rep(NA,length(y))
  y_impute_val = orig_data_vec
  
  if(length(test_inds)==0){
    return(cbind(y,y_impute_val))
  }
  
  n_test = length(test_inds)
  n_train = length(train_inds)
  d = distmat[test_inds,train_inds,drop=FALSE]
  NN = t(apply(d, 1, order))
  
  pred = apply(NN[, 1:k, drop=FALSE], 1, function(nn){
    tmp = rle(sort(y[train_inds[nn]]))  # rle(sort(c(3,4,4,2,5,7,3,3)))
    return(tmp$values[which.max(tmp$lengths)])
  })
  
  y[test_inds] = pred
  
  K = length(y_mus)
  unif_inds = which(pred==(K+1))
  
  if(length(unif_inds)>0){
    oth_inds = which(pred<=K)
    y_impute_val[test_inds[unif_inds]] = runif(length(unif_inds), 0, max_gene_expr)
    y_impute_val[test_inds[oth_inds]] =  y_mus[pred[oth_inds]]
  }else{
    y_impute_val[test_inds] = y_mus[pred]
  }
  
  return(cbind(y,y_impute_val))
}

# impute missing values (indicated by orig_data_vec==0) by kNN, may fail if all kNN are dropout
# NN: n x n matrix, n cells, nearest neighbours for each cell
# orig_data_vec: n x 1 vector, original (transformed) expression values for this gene
# gmm_label_vec: n x 1 vector, labels (component id) for this gene, 0 is dropout, others are 1 to K+1
# gmm_mus: K dimentional vector, each element is a component mean, dropout component (first) is removed, the last of which is the largest expression value for uniform component 
# k: k nearest neighbors
# return: y with all dropout labels (0) imputed, y_impute_val is imputed missing values (original values are copied)
gmm_impute1 = function(NN, orig_data_vec, gmm_label_vec, gmm_mus, k, test_inds=NA){
  y = gmm_label_vec
  y_mus = gmm_mus
  # n_cell = length(gmm_label_vec)
  max_gene_expr = max(orig_data_vec)
  
  if(is.na(test_inds)){
    test_inds = which(y==0)
  }else{
    stopifnot(all(gmm_label_vec[test_inds]==0))  # test_inds must not exceed range and must be dropout
  }

  y_impute_val = orig_data_vec
  
  if(length(test_inds)==0){
    return(cbind(y,y_impute_val))
  }
  
  pred = apply(NN[test_inds, 1:k+1, drop=FALSE], 1, function(nn){
    tmpinds = which(y[nn]!=0)
    if(length(tmpinds)==0){
      return(0)
    }else{
      nn = nn[tmpinds]
    }
    tmp = rle(sort(y[nn]))  # rle(sort(c(3,4,4,2,5,7,3,3)))
    return(tmp$values[which.max(tmp$lengths)])
  })
  
  test_inds = test_inds[pred!=0]
  if(length(test_inds)==0){
    return(cbind(y,y_impute_val))
  }
  
  y[test_inds] = pred
  
  K = length(y_mus)
  unif_inds = which(pred==(K+1))
  
  if(length(unif_inds)>0){
    oth_inds = which(pred<=K)
    y_impute_val[test_inds[unif_inds]] = runif(length(unif_inds), 0, max_gene_expr)
    y_impute_val[test_inds[oth_inds]] =  y_mus[pred[oth_inds]]
  }else{
    y_impute_val[test_inds] = y_mus[pred]
  }
  
  return(cbind(y,y_impute_val))
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

####  NOTE: specify the GMM result and gene express matrix here  #######
GMM_res = GMM_res_100  # a list, 1 x n_gene
gene_expr_mat = matrix(runif(315*100),315,100)   # n_cell x n_gene

data = gen_data(GMM_res, gene_expr_mat)
n_cell = data$n_cell
n_gene = data$n_gene

# filter genes
keep_rate_cutoff = 0.8
keep_rate = apply(data$X>0,2,sum)/n_cell
gene_inds = which(keep_rate>keep_rate_cutoff)

# perform imputation
Z = data$Z[,gene_inds]
mus = data$mus[gene_inds]
orig_data = data$X[,gene_inds]

# original imputation method
res = gmm_impute(orig_data, Z, mus)
X_impute = res$impute_label_mat
X_val_impute = res$impute_val_mat # both original and imputed values are there

# new method, may not impute if nearest neighbours all dropout
res1 = gmm_impute_v1(orig_data, Z, mus)
X_impute1 = res1$impute_label_mat
X_val_impute1 = res1$impute_val_mat # both original and imputed values are there
