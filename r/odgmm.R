### functions  ####
odgmm = function(X, # X: data, a 1-dimensional vector, must be log transformed of count vector, i.e. log(count+1)
                 n_max_comp=5, # maximum number of components
                 n_min_comp=3, # minimum number of components
                 mu0 = 0, # mean of dropout component
                 sigma0 = 0.1, # std of dropout component
                 min_variance = 1e-6,  # minimum variance allowed
                 min_beta = 0.9,  # mimimum beta allowed
                 min_ws = 0.01, # minimum Gaussian component weight
                 unif_max = NA, # maximum value for uniform component
                 debug=FALSE){
  
  # perform overdispersed Gaussian mixture modeling on input data X, with maximum K components
  
  # functions for transforming parameters
  alpha2mu = function(alpha){
    # global var mu0
    return ( cumsum(c(mu0,alpha)) )
  }
  mu2alpha = function(mus){
    return ( diff(mus) )
  }
  beta2sigma = function(beta){
    # beta is for the variance
    # global var sigma0
    return( cumprod( c(sigma0,sqrt(beta)) ) )
  }
  beta2var = function(beta){
    # global var sigma0
    return( cumprod(c(var0,beta)) )
  }
  sigma2beta = function(sigmas){
    return (exp(diff(log(sigmas)))^2)
  }
  var2beta = function(vars){
    return (exp(diff(log(vars))))
  }
  
  # functions for EM algorithms
  estep = function(X, alpha, beta, ws){
    mus = alpha2mu(alpha)
    sds = beta2sigma(beta)
    K = length(alpha)
    N = length(X)
    Z = matrix(0,nrow=N,ncol=K+2)
    XX = matrix(rep(X,K+1), nrow=N, ncol = K+1)
    XX = t(dnorm(t(XX), mean = mus, sd = sds, log = TRUE))
    Z[,seq(K+1)] = XX + matrix( rep( log(ws[seq(K+1)]), each=N ), nrow=N ) 
    Z[,K+2] = log(ws[K+2]) + uni_log_density
    
    if(any(is.na(Z))){
      print('Z contains NA')
    }
    
    maxz = apply(Z,1,max)
    Z = exp(Z - maxz)
    return (Z/apply(Z,1,sum))
  }
  cal_alpha = function(X,Z,alpha,beta,j){
    mus = alpha2mu(alpha)
    vars = beta2var(beta)
    
    # variance of some component might be 0, limit it to be at least min_variance
    vars[vars<min_variance] = min_variance
    
    K = length(alpha)
    N = length(X)
    XX = rep(0,N)
    ZZ = rep(0,N)
    for(k in seq(j,K))
    {
      tmp = Z[,(k+1)]/vars[k+1]
      XX = XX + (X - mus[k+1] + alpha[j])*tmp
      ZZ = ZZ + tmp
    }
    tsum = sum(ZZ)  # vars can be 0
    if(is.na(tsum) || tsum==0){
      return(NA)
    }else{
      # return (sum(XX)/tsum)
      ret_val = sum(XX)/tsum
      if(ret_val<min_alpha){
        return(min_alpha)
      }else{
        return(ret_val)
      }
      
    }
  }
  cal_beta = function(X,Z,alpha,beta,j){
    mus = alpha2mu(alpha)
    vars = beta2var(beta)
    
    # variance of some component might be 0, limit it to be at least min_variance
    vars[vars<min_variance] = min_variance
    
    K = length(alpha)
    N = length(X)
    XX = rep(0,N)
    for(k in seq(j,K))
    {
      XX = XX + Z[,(k+1)] * (X - mus[k+1])^2 * beta[j] / vars[k+1]
    }
    tsum = sum(Z[,seq(j+1,K+1)])
    if(is.na(tsum) || tsum==0){
      return(NA)
    }else{
      ret_val = sum(XX)/tsum
      if(ret_val<min_beta){
        return(min_beta)
      }else{
        return(ret_val)
      }
    }
  }
  cal_ws = function(Z){
    N = dim(Z)[1]
    ws = apply(Z,2,sum)/N
    return (ws)
  }
  mstep = function(X,Z,alpha,beta,ws){
    K = length(alpha)
    new_ws = cal_ws(Z)
    new_alpha = alpha
    new_beta = beta
    for(j in seq(K,1,-1))
    {
      new_alpha[j] = cal_alpha(X,Z,alpha,beta,j)
      new_beta[j] = cal_beta(X,Z,alpha,beta,j)
    }
    return(list(alpha=new_alpha,beta=new_beta,ws=new_ws))
  }
  exp_log_like = function(X, Z, alpha, beta, ws){
    mus = alpha2mu(alpha)
    sds = beta2sigma(beta)
    K = length(alpha)
    N = length(X)
    
    XX = matrix(rep(X,K+1), nrow=N, ncol = K+1)
    
    tmpX = t(dnorm(t(XX), mean = mus, sd = sds, log = TRUE))  # may be Inf if sd=0 (single data point)
    tmpX = tmpX + matrix( rep( log(ws[seq(K+1)]), each=N ), nrow=N) 
    tmpX = cbind(tmpX, rep(log(ws[K+2])+uni_log_density, N)) # log(0)=-Inf
    tmpX[Z==0] = 0  # handle case log(0)=-Inf in tmpX
    return(sum(tmpX*Z))
  }
  elbo = function(X,Z,alpha,beta, ws){
    LZ = Z
    LZ[Z!=0] = log(Z[Z!=0])
    entropy = -1 * Z * LZ
    
    lb = exp_log_like(X,Z,alpha,beta, ws)+sum(entropy)
    # if(is.na(lb) || lb==-Inf){
    #   # exp_log_like(X,Z,alpha,beta, ws)
    #   print("lower bounder is na or -Inf.")
    # }
    return (lb)
  }
  cal_bic = function(X,Z,alpha,beta,ws){
    N = length(X)
    K = length(alpha)
    res = -2*exp_log_like(X,Z,alpha,beta, ws) + (3*K+1)*log(N) # the smaller bic, the better model
    return(res)
  }
  # generate the result
  gen_res = function(alpha=NA,beta=NA,ws=NA,lb=NA,Z=NA,bic=NA,label=NA){
    if(any(is.na(alpha))){
      mus = NA
    }else{
      mus = alpha2mu(alpha)
    }
    if(any(is.na(beta))){
      sigmas = NA
    }else{
      sigmas = beta2sigma(beta)
    }
    # unif_max from outer function
    if(any(is.na(ws))){
      return(list(alpha=alpha,beta=beta,mus=mus,sigmas=sigmas,ws=ws,unif_ws=NA,unif_max=unif_max,lb=lb,bic=bic,label=label))
    }else{
      ncomp = length(ws)
      if(is.matrix(Z) && is.na(label)){
        label = apply(Z, 1, which.max)
        label[label==ncomp] = 0 # Uniform component is labeled as 0
      }
      return(list(alpha=alpha,beta=beta,mus=mus,sigmas=sigmas,ws=ws[-ncomp],unif_ws=ws[ncomp],unif_max=unif_max,lb=lb,bic=bic,label=label))
    }
  }
  disp_paras = function(res, infostr=''){
    cat(paste( infostr,
               paste0('alpha=',sprintf("%.3f", res$alpha),collapse=" "),
               paste0('beta=',sprintf("%.3f", res$beta),collapse=" "),
               paste0('ws=',sprintf("%.3f", res$ws),collapse=" "),
               paste0('mus=',sprintf("%.3f", alpha2mu(res$alpha)),collapse=" "),
               paste0('vars=',sprintf("%.3f", beta2var(res$beta)),collapse=" "),
               paste0('unif_ws=',sprintf("%.3f",  res$unif_ws, collapse=" ") ),
               paste0('unif_max=',sprintf("%.3f", res$unif_max, collapse=" ") ),
               paste0('lb=', sprintf("%.3f", res$lb),  collapse=""),
               paste0('bic=',sprintf("%.3f", res$bic),  collapse=""),
               '','',
               sep="\n" 
            ) 
        )
  }

  # perform em algorithm for given input parameters
  em_algo = function(X, init_mus, init_sgs, init_ws, nround=200, debug=FALSE){
    alpha = mu2alpha(init_mus)
    beta = var2beta(init_sgs)
    ws = init_ws
    Z = estep(X, alpha, beta, ws)
    lb = -Inf
    lb_arr = rep(NA,nround)
    curr_res = gen_res()
    for (i in seq(nround)){
      if(debug){
        print(paste0('iteration=',i,'  lb=',lb))
      }
      Z = estep(X, alpha, beta, ws)
      res = mstep(X,Z,alpha,beta,ws)
      
      if(any(is.na(res$alpha)) || any(is.na(res$beta)) || any(is.na(ws))){
        return( curr_res )
      }
      
      lb_new = elbo(X, Z, res$alpha, res$beta, res$ws)
      lb_arr[i] = lb_new
      alpha = res$alpha
      beta = res$beta
      ws = res$ws
      
      if(lb_new==-Inf || is.na(lb_new) || is.na(lb)){
        if(debug){
          cat(paste0('lb_new or lb is NaN. Use result in iteration ',i,'.\n'))
          cat("lb_new=",lb_new,"\n")
          cat("lb=",lb,"\n")
        }
        return(curr_res)
        # stop(" of lb_new or lb is NaN.")
      }
      
      if( abs(lb_new-lb) < abs(1e-6*lb) )
      {
        lb = lb_new
        break
      }
      lb = lb_new
      bic = cal_bic(X,Z,alpha,beta,ws)
      curr_res = gen_res(alpha,beta,ws,lb,Z,bic=bic)
    }
    if(i==nround){
      if(debug){
        print(paste0('Run all ',i,' iterations.',' lb=',lb))
      }
    }else{
      if(debug){
        print(paste0('Converge in ',i,' iterations.',' lb=',lb))
      }
    }
    # lb_arr = lb_arr[!is.na(lb_arr)]
    # if(debug){
    #   idx = seq(3,nround)
    #   plot(idx,lb_arr[idx])
    # }
    return(curr_res)
  }
  
  # run multiple runs and select the best one
  # X: raw data
  # K: number of Gaussian component, excluding the dropout component
  em_optim = function(X,K,debug=FALSE){
    stopifnot(length(X)>0)
    if(length(X)==1 && K==1){
      alpha = X
      beta = min_beta
      Z = matrix(c(0,1,0), nrow=1)
      ws=c(0,1,0)
      lb = elbo(X, Z, alpha, beta, ws)
      bic_val = cal_bic(X,Z,alpha,beta,ws)
      return(gen_res(alpha = alpha, 
                     beta = beta,
                     ws=ws, 
                     Z=Z,
                     lb = lb,
                     bic = bic_val
      ))
    }else if(length(X)==1){
      return( gen_res() )
    }
    
    # depends on global mu0, sigma0
    if(var(X)<=var0 && K>1){
      return( gen_res() )
    }else if(K==0 || var(X)<=var0){  # K=1, assign all X to this component
      xqu = mean(X)
      squ = c(var0, var(X))
      n_trial = 1
    }else{
      xqu = unname(quantile(X, probs = seq(0.05, 1, 0.05)))
      squ = seq(log(var0),log(var(X)),length.out=21)
      n_trial = 20
    }
    
    lb_arr = rep(-Inf,n_trial)
    bic_arr = rep(Inf,n_trial)
    res_list = vector("list",n_trial)
    
    for(i in seq(n_trial)){
      if(debug){
        cat('\n-----------------n_comp=',K+2,' | ', 'i_trial=',i, ' | n_trial=',n_trial,' -------------\n',sep='')
      }
      res_list[[i]]=tryCatch(
        {
          # initialization
          tmpw = runif(K+2)+1
          tmpw[1] = 0.5
          tmpw[length(tmpw)]=0.5
          init_ws = tmpw/sum(tmpw)
          init_mus = c(mu0, sort( sample(xqu,K)  ))
          init_sgs = c(var0, exp( sample(squ[-1],K,replace = TRUE))  )
          em_algo(X, init_mus, init_sgs, init_ws, debug=debug)
        },error=function(e){
          next()
        }
      )
      lb_arr[i] = res_list[[i]]$lb
      bic_arr[i] = res_list[[i]]$bic
    }
    min_ind = which.min(bic_arr)
    res = res_list[[min_ind]]
    
    cat('\n-----------------n_comp=',K+2,' final result -------------\n',sep='')
    if(!is.null(res)){
      if(K==0){
        print('K=0')
      }
      # display result
      disp_paras(res, infostr='')
    }else{
      res = gen_res()
      cat("No results available for K=",K+2,"\n\n")
    }
    
    return(res)
  }
  
  # main code for odgmm
  # initialization
  # mu0,sigma0
  var0 = sigma0^2
  min_alpha = 3*sigma0 # minimum alpha allowed
  
  stopifnot(n_min_comp>=3)
  
  # assign all 0 directly to the first component
  n_data_point = length(X)
  non_zero_inds = which(X!=0)
  zero_inds = which(X==0)
  n_zero_point = length(zero_inds)
  non_zero_weight = length(non_zero_inds)/n_data_point
  zero_weight = 1-non_zero_weight
  
  # output information here
  infostr = sprintf('%d/%d=%.1f%% are zeros, directly assigned to dropout component.\n',
                    n_zero_point, n_data_point, n_zero_point/n_data_point*100)
  cat(infostr,'\n')
  
  # in case input data are all zeros
  if(length(non_zero_inds)==0){
    cat('all input data are zero. Quit.\n')
    lab_vec = rep(1,n_data_point)
    return(gen_res(alpha=min_alpha, beta = min_beta, ws=c(1,0), label=lab_vec, unif_max=unif_max, unif_ws=0))
  }
  
  X = X[non_zero_inds]
  
  stopifnot(n_max_comp>=n_min_comp)
  if(n_min_comp<2){
    stop(paste0("n_min_comp=",n_min_comp," must at least be 2. One dropout component and one signal component."))
  }

  if(is.na(unif_max)){
    unif_max = ceiling(max(X))
  }
  uni_log_density = -log(unif_max)
  
  lb_arr = rep(NA,n_max_comp)
  bic_arr = rep(NA,n_max_comp)
  res_list = vector("list",n_max_comp)
  flag_arr = rep(F,n_max_comp)
  ws_flag_arr = rep(F,n_max_comp)
  
  for(i in seq(n_max_comp,n_min_comp)){
    K = i - 1
    cat('\n********************* n_component=',i, ' **********************\n',sep='')
    res_list[[i]] = em_optim(X, K-1,debug=debug)
    lb_arr[i] = res_list[[i]]$lb
    bic_arr[i] = res_list[[i]]$bic
    
    if(!is.na(res_list[[i]]$lb)){
      if(K>1 && all(res_list[[i]]$ws[2:K]>=min_ws)){
        ws_flag_arr[i] = TRUE
      }
      if(K==0){
        ws_flag_arr[i] = TRUE
      }
      flag_arr[i] = TRUE
    }
  }
  if(n_max_comp>2 && any(ws_flag_arr[3:n_max_comp])){
    valid_inds = which(ws_flag_arr)
  }else{
    valid_inds = which(flag_arr)
  }
  
  if(length(valid_inds)>0){
    min_ind = which.min(bic_arr[valid_inds])
    min_ind = valid_inds[min_ind]
    res = res_list[[min_ind]]
  }else{
    stop('No valid results available.\n')
  }
  
  label = rep(1,n_data_point)
  label[non_zero_inds] = res$label
  res$label =label
  k = length(res$ws)
  tmp_ws = c(res$ws,res$unif_ws)
  tmp_ws = non_zero_weight * tmp_ws
  tmp_ws[1] = tmp_ws[1] + zero_weight
  tmp_ws = tmp_ws/sum(tmp_ws)
  res$ws = tmp_ws[1:k]
  res$unif_ws = tmp_ws[k+1]
  
  cat('\n****************** Final Result ',k,' Gaussian components ****************\n',sep='')
  infostr = sprintf('%d/%d=%.1f%% zeros are added to and renormalized in the final result.',
                    n_zero_point, n_data_point, n_zero_point/n_data_point*100)
  disp_paras(res, infostr=infostr)
  
  return(res)
}

gmmpdf = function(x, mus, sigmas, ws, log=FALSE){
  K = length(mus)
  N = length(x)
  y = rep(0,N)
  for(k in seq(K)){
    y = y + ws[k] * dnorm(x,mus[k],sigmas[k])
  }
  if(log){
    return(log(y))
  }else{
    return (y)
  }
}

odgmm_diff = function(res1, res2){
  # remove dropout component, extend uniform component to given range
  # res is odgmm output
  renorm = function(res, unif_max){
    n = length(res$ws)
    unif_ws = ifelse(res$unif_ws>1e-3,res$unif_ws,1e-3)
    ws = res$ws[2:n]
    
    tmpws = sum(c(unif_ws,ws))
    unif_ws = unif_ws/tmpws
    ws = ws/tmpws
    
    mus = res$mus[2:n]
    sigmas = res$sigmas[2:n]
    return(list(ws=ws,mus=mus,sigmas=sigmas,unif_max=unif_max,unif_ws=unif_ws))
  }
  sample_normal = function(n, mu, sigma, min_x, max_x){
    x = rnorm(n, mu, sigma)
    valid_inds = x>=min_x & x<=max_x
    return(x[valid_inds])
  }
  # sample from gmm+uniform mixture distribution
  # should not generate value lower than 0
  sample_gmm_unif_mix = function(model, n){
    m = length(model$ws)
    reslist = vector("list",m+1)
    for(i in seq(m)){
      tmpn = round(n*model$ws[i])
      reslist[[i]] = sample_normal(tmpn, model$mus[i], model$sigmas[i], 0, model$unif_max)
    }
    tmpn = round(n*model$unif_ws)
    reslist[[m+1]] = rnorm(tmpn)*model$unif_max
    return(unlist(reslist))
  }
  mix_pdf = function(x, model){
    res = gmmpdf(x, model$mus, model$sigmas, model$ws)
    res = res+model$unif_ws/model$unif_max
    return(res)
  }
  kl_div = function(x, model1, model2){
    y1 = log(mix_pdf(x,model1))
    y2 = log(mix_pdf(x,model2))
    res = sum(y1-y2)/length(x)
    return(res)
  }
  max_x = max(res1$unif_max, res2$unif_max)
  model1 = renorm(res1, max_x)
  model2 = renorm(res2, max_x)
  
  nsample = 2000
  x1 = sample_gmm_unif_mix(model1, nsample)
  x2 = sample_gmm_unif_mix(model2, nsample)
  
  dist12 = kl_div(x1, model1, model2)
  dist21 = kl_div(x2, model2, model1)
  
  return(list(dist12=dist12, dist21=dist21, model1=model1, model2=model2))
}

# make a plot for gmm+unif mixture model
plot_mix_model = function(model, col='red', lty=1, lwd=2){
  emus = model$mus
  esgs = model$sigmas
  ews = model$ws
  
  x = seq(-2,model$unif_max+1,0.01)
  ey = gmmpdf(x, emus, esgs, ews) + model$unif_ws/model$unif_max
  lines(x,ey,col=col,lty=lty, lwd=lwd)
}

require(rhierbaps)
# Z is n_cell x n_gene matrix, all elements range 0~4
# 0 stands for uniform component, 1 stands for dropout component
# returns a data.frame which contains the partition of n_max_depth=2 layers
odgmm_cluster = function(Z, n_max_cluster, n_max_depth=2, cell_names=NULL){
  n_cell = dim(Z)[1]
  n_gene = dim(Z)[2]
  snp_mat = as.factor(Z)
  levels(snp_mat) = c('-','a','c','g','t')
  snp_mat = matrix(snp_mat, nrow=n_cell, ncol = n_gene)
  
  if(is.null(cell_names)){
    rownames(snp_mat) = seq(n_cell)
  }else{
    rownames(snp_mat) = cell_names
  }
  res = hierBAPS(snp_mat, max.depth = n_max_depth, n.pops = n_max_cluster, quiet = TRUE)
  return(res$partition.df[,2:(n_max_depth+1)])
}

# imputation
# cluster_label_vec: cluster labels for each cell, n_cell x 1 vector
# X: log transformed original data, n_cell x n_gene
# Z: odgmm component label for each gene, n_cell x n_gene
# mus_list: mus of each gene, n_gene x 1 list
# cutoff: impute a gene if the proportion of cells that has values (non-dropout or uniform) exceeds this cutoff
# output an imputed matrix
odgmm_impute = function(cluster_label_vec, X, Z, mus_list, cutoff=0.8, debug=FALSE){
  impute_vec = function(label_vec, value_vec, mus){
    impute_val_vec = value_vec
    test_inds = which(label_vec==1)   # dropout component, 0 is uniform component
    train_inds = which(label_vec>1)   # real signals
    if(length(test_inds)==0) return(value_vec)
    pred_lab = label_vec[sample(train_inds,length(test_inds),replace = TRUE)]
    impute_val_vec[test_inds] = mus[pred_lab]
    return(impute_val_vec)
  }

  I = Z>1   # 0-uniform component, 1-dropout component
  uniq_lab= sort(unique(cluster_label_vec))
  n_pop = length(uniq_lab)
  n_gene = dim(X)[2]
  
  for(i in uniq_lab){
    tmpinds = which(cluster_label_vec==i)
    tmpn = length(tmpinds)
    tmp_p_sum = apply(I[tmpinds,,drop=FALSE],2,sum)/tmpn
    
    if(tmpn<2) next()
    
    tmp_impute_gene_inds = which(tmp_p_sum>=cutoff)
    for(g in tmp_impute_gene_inds){
      if(debug){
        cat(sprintf("impute cluster i=%d gene g=%d",i,g),"\n")
      }
      X[tmpinds, g] = impute_vec(Z[tmpinds,g], X[tmpinds,g], mus_list[[g]])
    }
  }
  return(X)
}

check_diff_res = function(ind, dist12_arr, dist21_arr, diff_res_list, gene_exp_log_mat1, gene_exp_log_mat2){
  tmpres = diff_res_list[[ind]]
  m1 = tmpres$model1
  m2 = tmpres$model2
  max_y = max(c(density(gene_exp_log_mat1[,ind])$y,
                density(gene_exp_log_mat2[,ind])$y ))*2
  plot(density(gene_exp_log_mat1[,ind],bw='sj'),ylim=c(0,max_y),
       main=paste0("dist12=",round(dist12_arr[ind],digits = 2)," dis21=",round(dist21_arr[ind],digits = 2)))
  lines(density(gene_exp_log_mat2[,ind],bw='sj'),lty=2)
  plot_mix_model(m1)
  plot_mix_model(m2,col="blue",lty=2)
}

####  main code #####
set.seed(1)

# ## part 1: simulate data from Gaussian mixture and perform ODGMM
# N = 2000
# wghs = c(0.2,0.2,0.6)
# components = sample(1:3,prob=wghs,size=N,replace=TRUE)
# mus = c(0, 2,  5)
# sds = c(0.1, 0.5, 1)
# sgs = sds^2
# 
# X1 = rnorm(n=N,mean=mus[components],sd=sds[components])
# res1 = odgmm(X1,5,debug=F) # only need to keep this result!!!
# # display fitting with original density
# plot(density(X1,bw='sj'),ylim=c(0,0.9),main=paste0("GMM vs kernel density"))
# plot_mix_model(res1)
# 
# N = 2000
# wghs = c(0.2,0.6,0.2)
# components = sample(1:3,prob=wghs,size=N,replace=TRUE)
# mus = c(0, 2,  5)
# sds = c(0.1, 0.5, 1)
# sgs = sds^2
# 
# X2 = rnorm(n=N,mean=mus[components],sd=sds[components])
# res2 = odgmm(X2,5,debug=F) # only need to keep this result!!!
# # display fitting with original density
# plot(density(X2,bw='sj'),ylim=c(0,0.9),main=paste0("GMM vs kernel density"))
# plot_mix_model(res2)
# 
# ## part 2: perform differential gene expression
# # example 2.1
# # the KL divergence is not symmetric, so there are two values
# tmpres = odgmm_diff(res1,res2)
# m1 = tmpres$model1
# m2 = tmpres$model2
# 
# cat('KL divergence of m1 to m2 dist12=',tmpres$dist12,"\n")
# cat('KL divergence of m2 to m1 dist21=',tmpres$dist21,"\n")
# 
# plot(density(X1,bw='sj'),ylim=c(0,0.9),main=paste0("Kernel density and gmm+unif fits"))
# plot_mix_model(m1)
# lines(density(X2,bw='sj'),lty=2)
# plot_mix_model(m2,col="blue",lty=2)
# # cannot get legend properly plotted
# # legend("topleft",
# #        c("X1 kernel density", "X1 gmm+unif fit, exluding dropout",
# #                             "X2 kernel density", "X2 gmm+unif fit, exluding dropout"),
# #        adj = c(0, 0.1),
# #        )
# #        # col=c("black","red", "black","blue"),
# #        # lty=c(1,1,2,2))
# 
# # example 2.2, perform differential expression to select genes
# 
# gene_exp_mat = t(readRDS("expr_100genes_1000cells.rds"))  # n_cell x n_gene
# gene_exp_mat = as.matrix(gene_exp_mat)
# n_cell = dim(gene_exp_mat)[1]
# n_gene = dim(gene_exp_mat)[2]
# gene_exp_log_mat = log(gene_exp_mat+1)
# 
# # generate artifical groups
# group_vec = sample(2,n_cell,replace = T, prob = c(0.5,0.5))
# 
# # filter gene by expression proportion
# min_n_cell = 10 # at least min_n_cell cells with expression
# valid_gene_inds = which(apply(gene_exp_log_mat>0, 2, sum)>min_n_cell)
# gene_exp_log_mat = gene_exp_log_mat[,valid_gene_inds]
# n_gene = length(valid_gene_inds)
# 
# # separate into two groups
# gene_exp_log_mat1 = gene_exp_log_mat[group_vec==1,]
# gene_exp_log_mat2 = gene_exp_log_mat[group_vec==2,]
# res_list1 = vector("list",n_gene)
# res_list2 = vector("list",n_gene)
# diff_res_list = vector("list",n_gene)
# dist12_arr = rep(NA,n_gene)
# dist21_arr = rep(NA,n_gene)
# 
# keep_rate1 = apply(gene_exp_log_mat1>0, 2, sum)/n_cell
# keep_rate2 = apply(gene_exp_log_mat2>0, 2, sum)/n_cell
# 
# # odgmm for each gene, and differential expression
# for(i in seq(1,n_gene)){
#   cat('odgmm for gene i=',i,"\n")
#   tmpmax = ceiling( max( c(gene_exp_log_mat1[,i],gene_exp_log_mat2[,i]) ) )
#   res_list1[[i]] = odgmm(gene_exp_log_mat1[,i],5,unif_max=tmpmax, debug=F)
#   res_list2[[i]] = odgmm(gene_exp_log_mat2[,i],5,unif_max=tmpmax, debug=F)
#   diff_res_list[[i]] = odgmm_diff(res_list1[[i]],res_list2[[i]])
#   dist12_arr[i] = diff_res_list[[i]]$dist12
#   dist21_arr[i] = diff_res_list[[i]]$dist21
# }
# 
# # rank the distance the choose the significant ones
# min_sig_rate = 0.3
# keep_rate_sig_gene_inds = which(abs(keep_rate1 - keep_rate2)>min_sig_rate)
# min_keep_rate = 0.1 
# # only for genes without significant dropout difference and high enough keep rate
# diff_valid_gene_inds = which(abs(keep_rate1 - keep_rate2)<=min_sig_rate &
#                          keep_rate1>min_keep_rate & keep_rate2>min_keep_rate)
# diff_dist_arr = c(dist12_arr[diff_valid_gene_inds], dist21_arr[diff_valid_gene_inds])
# top_prop = 0.05
# diff_dist_cutoff = quantile(diff_dist_arr,1-top_prop)
# diff_sig_gene_inds = diff_valid_gene_inds[ dist12_arr[diff_valid_gene_inds]>=diff_dist_cutoff
#                                            | dist21_arr[diff_valid_gene_inds]>=diff_dist_cutoff ]
# 
# keep_diff_sig_gene_inds = which(abs(keep_rate1 - keep_rate2)>min_sig_rate &
#                             (dist12_arr>diff_dist_cutoff | dist21_arr>diff_dist_cutoff) )
# 
# cat("Significant KL divergence (top 5%, similar keep rate) gene inds: ",diff_sig_gene_inds, '\n')
# cat("Significant different keep rate gene inds: ",keep_rate_sig_gene_inds, '\n')
# cat("Significant different keep rate & KL divergence gene inds: ",keep_diff_sig_gene_inds, '\n')
# 
# # plot the significant genes
# check_diff_res(diff_sig_gene_inds[1], dist12_arr,dist21_arr, diff_res_list, gene_exp_log_mat1, gene_exp_log_mat2)

# ## part 3: clustering gene expression matrix
# # example 3.1
# # res is a n_cell x 2 matrix, 1 column is the first layer clustering, 2nd column is second layer clustering
# Z = matrix(round(runif(600)*4), nrow=60, ncol=10)
# res = odgmm_cluster(Z,20)
# 
# # example 3.2
# gene_exp_mat = t(readRDS("expr_100genes_1000cells.rds"))  # n_cell x n_gene
# gene_exp_mat = as.matrix(gene_exp_mat)
# n_cell = dim(gene_exp_mat)[1]
# n_gene = dim(gene_exp_mat)[2]
# gene_exp_log_mat = log(gene_exp_mat+1)
# 
# res_list = vector("list",n_gene)
# Z = matrix(0, nrow = n_cell, ncol = n_gene)
# 
# # perform odgmm for each gene
# for(i in seq(1,10)){
#   cat('odgmm for gene i=',i,"\n")
#   res_list[[i]] = odgmm(gene_exp_log_mat[,i],5, debug=F)
#   Z[,i] = res_list[[i]]$label
# }
# 
# # filter genes that has low counts
# min_n_cell = 5
# valid_gene_inds = which(apply(gene_exp_log_mat>0, 2, sum)>min_n_cell)
# ZZ = Z[,valid_gene_inds]
# 
# # perform clustering
# partition = odgmm_cluster(ZZ,20)


# ## part 4: perform imputation
# gene_exp_mat = t(readRDS("expr_100genes_1000cells.rds"))  # n_cell x n_gene
# gene_exp_mat = as.matrix(gene_exp_mat)
# n_cell = dim(gene_exp_mat)[1]
# n_gene = dim(gene_exp_mat)[2]
# 
# # log transformation to the data
# gene_exp_log_mat = log(gene_exp_mat+1)
# 
# # perform odgmm for each gene 
# res_list = vector("list",n_gene)
# for(i in seq(1,n_gene)){
#   cat('odgmm for gene i=',i,"\n")
#   res_list[[i]] = odgmm(gene_exp_log_mat[,i],5,debug=F)
# }
# 
# # perform clustering to the data
# Z = matrix(0, nrow=n_cell, ncol=n_gene)
# mus_list = vector("list", n_gene)
# for(i in seq(n_gene)){
#   Z[,i] = res_list[[i]]$label
#   mus_list[[i]] = res_list[[i]]$mus
# }
# partition = odgmm_cluster(Z,20)
# 
# # perform imputation
# imp_data_mat = odgmm_impute(partition[,2], gene_exp_log_mat, Z, mus_list, cutoff=0.8, debug=T)

















