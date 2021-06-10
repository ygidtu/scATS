# ats mixture model inference
atsmix = function(n_max_ats, # maximum number of ATS sites
                  n_min_ats=1, # minimum number of ATS sites
                  
                  # information for all DNA fragments, must be n_frag x 1 vector if specified
                  st_arr, # l, start location of each DNA fragment, from 5' to 3'on the 5'-UTR
                  en_arr, # r, end location of each DNA fragment

                  L=2000,  # length of UTR region
                  
                  # fragment size information
                  mu_f=300, # fragment length mean
                  sigma_f=50, # fragment length standard deviation
                  
                  # pa site information
                  min_ws = 0.01, # minimum weight of ATS site
                  max_unif_ws = 0.1, # maximum weight of uniform component
                  max_beta = 50, # maximum std for ATS site
                  
                  # inference with fixed parameters
                  fixed_inference_flag = FALSE,
                  
                  #single end mode
                  single_end_mode = FALSE, 
                  ){
  
  cal_z_k = function(alpha_arr, beta_arr, ws, k, log_zmat){
    K = length(ws) - 1  # last component is uniform component
    if(k<=K){
      log_zmat[,k] = log(ws[k]) + lik_l_ab(st_arr, alpha_arr[k], beta_arr[k], log=T)
    }else{
      log_zmat[,K+1] = log(ws[K+1]) + unif_log_lik
    }
    return(log_zmat)
  }
  
  # Z is log likelihood
  norm_z = function(Z){
    tmp_max_vec = apply(Z,1,max)
    Z = exp( Z - tmp_max_vec )
    Z = Z/apply(Z,1,sum)
    
    return(Z)
  }

  # maximize ws given Z
  maximize_ws = function(Z){
    ws = apply(Z,2,sum)/dim(Z)[1]
    m = length(ws)
    
    if(ws[m]>max_unif_ws){
      ws[-m] = (1-max_unif_ws)*ws[-m]/sum(ws[-m])
      ws[m] = max_unif_ws
    }
    
    return( ws )
  }
  
  mstep = function(alpha_arr, beta_arr, ws, Z, k){
    
    K = length(ws) - 1 # last component is uniform component
    new_ws = maximize_ws(Z)
    new_alpha_arr = alpha_arr
    new_alpha_arr[k] = sum(Z[,k]*st_arr)/sum(Z[,k])
    
    new_beta_arr = beta_arr;
    tmp_beta = sqrt( sum(Z[,k]*(st_arr-new_alpha_arr[k])^2)/sum(Z[,k]) )
    if( is.na(tmp_beta) ){
      tmp_beta = max_beta
    }else if(tmp_beta<step_size){
      tmp_beta = step_size
    }
    # cat("predef_beta_arr=",predef_beta_arr,"\n")
    
    new_beta_arr[k] = predef_beta_arr[which.min(abs(predef_beta_arr-tmp_beta))]
    if(new_beta_arr[k]>max_beta){
      new_beta_arr[k] = max_beta
    }
    
    return(list(alpha_arr=new_alpha_arr, beta_arr=new_beta_arr, ws=new_ws))
  }
  
  elbo = function(log_zmat, Z){
    LZ = Z
    LZ[Z!=0] = log(Z[Z!=0])
    entropy = -1 * Z * LZ
    
    lb = exp_log_lik(log_zmat, Z) + sum(entropy)
    if(is.na(lb)){
      stop("lower bounder is na.")
    }
    return (lb)
  }
  
  # calculate the expected log joint likelihood
  exp_log_lik = function(log_zmat, Z){
    ZZ = Z * log_zmat
    ZZ[Z==0] = 0
    
    return(sum(ZZ))
  }
  
  # uniform component likelihood
  lik_f0 = function(log=F){
    if(log){
      return(-2*log(L))
    }else{
      return(1/L/L)
    }
  }
  
  # uniform component likelihood for single end case
  lik_f0_single = function(log=F){
    if(log){
      return( -log(L)+dnorm(mu_f,mean=mu_f,sd=sigma_f,log=log) )
    }else{
      return( 1/L*dnorm(mu_f,mean=mu_f,sd=sigma_f,log=log) )
    }
  }
  
  # p(l,r|alpha,beta)
  lik_lr_ab = function(l_arr, r_arr, alpha, beta, log=F){
    if(log){
      return(lik_l_ab(l_arr,alpha, beta, log=T)+lik_r_l(l_arr, r_arr, log=T))
    }else{
      return(ws * lik_l_ab(l_arr,alpha, beta) * lik_r_l(l_arr, r_arr))
    }
  } 
  
  # p(l|alpha, beta)
  lik_l_ab = function(l_arr, alpha, beta, log=F){
    return( dnorm(l_arr, mean=alpha, sd=beta, log=log)  )
  }
  
  # p(r|l)
  lik_r_l = function(l_arr, r_arr, log=F){
    return( dnorm(l_arr-r_arr, mean=mu_f, sd=sigma_f, log=log)  )
  }
  
  # generate random k such that each K is a group and no consecutive elements are the same
  gen_k_arr = function(K, n){
    if(K==0){
      return(rep(0,n))
    }
    if(K==1){
      return(rep(1,n))
    }
    if(K==2){
      nn = ceiling(n/K)
      res = rep(c(1,2),nn)
      return(res[seq(n)])
    }
    nn = ceiling(n/K)
    res = rep(0,nn*K)
    res[seq(1,K)]=sample(K)
    for(i in seq(2,nn)){
      st = (i-1)*K
      res[seq(st+1,st+K)]=sample(K)
      if(res[st]==res[st+1]){
        tmpind = sample(seq(2,K),1)
        tmp = res[st+tmpind]
        res[st+tmpind] = res[st+1]
        res[st+1] = tmp
      }
    }
    return(res[seq(n)])
  }
  
  cal_bic = function(log_zmat, Z){
    N = dim(Z)[1]
    K = dim(Z)[2]-1
    
    res = -2*exp_log_lik(log_zmat, Z) + (3*K+1) *log(N)  # the smaller bic, the better model
    return(res)
  }
  
  # perform inference given alpha_arr and beta_arr
  fixed_inference = function(alpha_arr, beta_arr){
    lb = -Inf
    lb_arr = rep(NA,nround)
    N = n_frag
    K = length(alpha_arr)
    
    ws = runif(K+1)
    ws[1:K] = ws[1:K]+1
    ws = ws/sum(ws)

    k_arr = gen_k_arr(length(alpha_arr), nround)
    
    log_zmat = matrix( 0, nrow=N, ncol=K+1 )
    for(k in seq(K+1)){
      log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k, log_zmat)
    }
    
    for(i in seq(nround)){

      # estep
      log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k_arr[i], log_zmat)
      Z = norm_z(log_zmat)
      
      # mstep, alpha and beta are fixed now
      ws = maximize_ws(Z)  
      
      lb_new = elbo(log_zmat, Z)
      lb_arr[i] = lb_new
      
      if(lb_new==-Inf){
        lb = -Inf
        break
      }
      
      if( abs(lb_new-lb) < abs(1e-6*lb) )
      {
        break
      }else{
        lb = lb_new
      }
    }

    bic = cal_bic(log_zmat, Z)
    label = apply(log_zmat, 1, which.max)
    
    lb_arr = lb_arr[!is.na(lb_arr)]

    return(list(ws=ws, alpha_arr=alpha_arr, beta_arr=beta_arr, lb_arr=lb_arr, bic=bic, label=label))
  }
  
  # perform inference for K components
  # theta_arr => theta_win_mat
  em_algo = function(ws, para_mat){
    lb = -Inf
    lb_arr = rep(NA,nround)
    N = n_frag
    K = length(ws)-1
    
    alpha_arr = para_mat[,1]
    beta_arr = para_mat[,2]
    
    k_arr = gen_k_arr(length(alpha_arr), nround)
    log_zmat = matrix( 0, nrow=N, ncol=K+1 )
    for(k in seq(K+1)){
      log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k, log_zmat)
    }
    
    for(i in seq(nround)){

      # estep
      log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k_arr[i], log_zmat)
      Z = norm_z(log_zmat)

      res = mstep(alpha_arr, beta_arr, ws, Z, k_arr[i])

      alpha_arr = res$alpha_arr
      beta_arr = res$beta_arr
      ws = res$ws
      
      lb_new = elbo(log_zmat, Z)
      lb_arr[i] = lb_new
      
      if(lb_new==-Inf){
        lb = -Inf
        break
      }
      
      if( abs(lb_new-lb) < abs(1e-6*lb) )
      {
        break
      }else{
        lb = lb_new
      }
    }
    
    bic = cal_bic(log_zmat, Z)
    
    lb_arr = lb_arr[!is.na(lb_arr)]

    sorted_inds = order(alpha_arr)
    alpha_arr = alpha_arr[sorted_inds]
    beta_arr = beta_arr[sorted_inds]
    ws[1:K] = ws[sorted_inds]
    return(list(ws=ws, alpha_arr=alpha_arr, beta_arr=beta_arr, lb_arr=lb_arr, bic=bic))
  }
  
  # calculate coverage 
  split_data = function(l_arr){
    get_boarder = function(sign_arr){
      tmp = rle(sign_arr)
      tmp_lens = tmp$lengths
      tmp_vals = tmp$values
      chng = cumsum( tmp_lens )
      tmpn = length(chng)
      st_arr = c(0,chng[-tmpn])+1
      en_arr = chng
      return(list(st_arr=st_arr, en_arr=en_arr))
    }
    
    # handle zero in sign array, half changed to previous/next segment non-zero sign 
    handle_zero = function(sign_arr){
      if(!any(sign_arr==0)){
        return(sign_arr)
      }
      
      # first sign is always 1
      stopifnot(sign_arr[1]==1)
      
      res = get_boarder(sign_arr)
      st_arr = res$st_arr
      en_arr = res$en_arr
      tmp_vals = rle(sign_arr)$values
      tmpn = length(tmp_vals)
      zero_inds = which(tmp_vals==0)
      for(i in zero_inds){
        st = st_arr[i]
        en = en_arr[i]
        if(i==tmpn | (en-st)<=1){
          sign_arr[st:en] = tmp_vals[i-1]
        }else{
          mid = round((st+en)/2)
          sign_arr[st:mid] = tmp_vals[i-1]
          sign_arr[(mid+1):en] = tmp_vals[i+1]
        }
      }
      return(sign_arr)
    }
    
    ks_res = density(x=l_arr,bw=5*step_size, kernel="gaussian",
                      n=1024) # step_size defined in outer function
    
    # change points
    sign_arr = sign(diff(c(-1,ks_res$y)))
    sign_arr = handle_zero(sign_arr)
    
    chng = cumsum( rle(sign_arr)$lengths )
    mode_arr = ks_res$x[ chng[seq(1,length(chng),2)] ]
    n_mode = length(mode_arr)
    boarder_arr = ks_res$x[ chng[seq(2,length(chng),2)] ][1:(n_mode-1)]
    st_arr = c(max(1,min(l_arr)-10),boarder_arr+1)
    en_arr = c(boarder_arr, max(l_arr)+10)
    ws_arr = rep(1,n_mode)
    for(i in seq(n_mode)){
      ws_arr[i] = sum(st_arr[i]<=l_arr & l_arr<=en_arr[i])
    }
    ws_arr = ws_arr/sum(ws_arr)
    
    return(list(n_win=n_mode, st_arr=st_arr, en_arr=en_arr, ws_arr=ws_arr, mode_arr=mode_arr))
  }
  
  init_para = function(n_ats){
    tmp_res = split_data(st_arr)
    if(tmp_res$n_win>=n_ats){
      tmpinds = sample(tmp_res$n_win, n_ats)
    }else{
      tmpinds = c(seq(tmp_res$n_win),sample(tmp_res$n_win, n_ats-tmp_res$n_win, replace = T, prob=tmp_res$ws))
    }
    tmpinds = sort(tmpinds)
    
    alpha_arr = rep(0, n_ats)
    for(i in seq(n_ats)){
      ti = tmpinds[i]
      tmpst = tmp_res$st_arr[ti]
      tmpen = tmp_res$en_arr[ti]
      alpha_arr[i] = runif(1, tmpst, tmpen)
    }

    beta_arr = runif(n_ats,step_size,max_beta)
    return(cbind(alpha_arr, beta_arr))
  }
  
  em_optim0 = function(n_ats){
    n_trial=20
    
    lb_arr = rep(-Inf,n_trial)
    bic_arr = rep(Inf,n_trial)
    res_list = vector("list",n_trial)
    
    for(i in seq(n_trial)){
      res_list[[i]] = tryCatch({
        ws = runif(n_ats+1)+1  # last component is for uniform component
        # ws[n_ats+1] = ws[n_ats+1]-1
        ws[n_ats+1] = max_unif_ws
        ws = ws/sum(ws)
        
        # initilize alpha_arr and beta_arr, considered pre_alpha_arr and pre_beta_arr
        para_mat = init_para(n_ats)
        
        em_algo(ws, para_mat)
      # res_list[[i]] = em_algo(ws, para_mat, debug=debug)
      },error=function(e){
        list(ws=NULL,alpha_arr=NULL,beta_arr=NULL,lb_arr=-Inf,bic=Inf)
      })
      if(is.null(res_list[[i]])){
        next()
      }
      lb_arr[i] = tail(res_list[[i]]$lb_arr,1)
      bic_arr[i] = res_list[[i]]$bic
    }
    
    min_ind = which.min(bic_arr)
    res = res_list[[min_ind]]
    
    return(res)
  }
  
  # remove components with weight less than min_ws
  rm_component = function(res, min_ws){
    K = length(res$alpha_arr)
    res$alpha_arr = round(res$alpha_arr)
    res$beta_arr = ceiling(res$beta_arr)
    rm_inds = which(res$ws[1:K]<min_ws)
    
    if(length(rm_inds)==0){
      return(res)
    }else{
      alpha_arr = res$alpha_arr[-rm_inds]
      beta_arr = res$beta_arr[-rm_inds]
       
      res = fixed_inference(alpha_arr,beta_arr)
      return(res)
    }
  }

  n_frag = length(st_arr)
  stopifnot(n_frag==length(en_arr))
  nround=50
  
  step_size=5
  if(max_beta<step_size){
    cat("max_beta=",max_beta," step_size=",step_size, "\n")
    stop('max_beta has to be greater than step_size.\n')
  }
  predef_beta_arr = seq(step_size,max_beta,step_size)

  lb_arr = rep(-Inf,n_max_ats)
  bic_arr = rep(Inf,n_max_ats)
  res_list = vector("list",n_max_ats)
  
  if(single_end_mode){
    st_arr = en_arr + mu_f
    unif_log_lik = lik_f0_single(log = T)
  }else{
    unif_log_lik = lik_f0(log = T)
  }
  
  for(i in seq(n_max_ats,n_min_ats)){
    res_list[[i]] = em_optim0(i)
    lb_arr[i] = tail(res_list[[i]]$lb_arr,1)
    bic_arr[i] = res_list[[i]]$bic
  }

  min_ind = which.min(bic_arr)
  res = res_list[[min_ind]]
  
  if(is.null(res$ws)){
    cat('\n Inference failed. No results available. \n')
    return(res)
  }
    
  # remove low weight component
  res = rm_component(res, min_ws)
  
  # calculate read assignment
  if(!'label' %in% names(res)){
    N = n_frag
    K = length(res$ws)-1
    log_zmat = matrix( 0, nrow=N, ncol=K+1 )
    for(k in seq(K+1)){
      log_zmat = cal_z_k(res$alpha_arr, res$beta_arr, res$ws, k, log_zmat)
    }
    # Z = norm_z(log_zmat)
    label = apply(log_zmat, 1, which.max)
    res$label = label
  }
  
  if(single_end_mode){
    res$beta_arr = res$beta_arr + sigma_f 
  }
  
  # check fragment size
  res$fragment_size_flag = 'normal'
  if(abs(mean(st_arr-en_arr+1)-mu_f)>2*sigma_f){
    warning(
      paste0('mean of fragment size is ',round(mean(st_arr-en_arr)),
             ', which deviates from mu_f=',mu_f, ' by more than 2*max_beta=', 2*sigma_f ,'bp.\n')
    )
    if(mean(st_arr-en_arr+1)-mu_f>0){
      res$fragment_size_flag = 'short'
      warning('fragment size shorter than expected.\n')
    }else{
      res$fragment_size_flag = 'long'
      warning('fragment size longer than expected. Potential exon skipping. \n')
    }
  }
  
  return(res)
}

est_density = function(alpha_arr, beta_arr, ws, L, x_arr){
  K = length(alpha_arr)
  y_arr = rep(0,length(x_arr))
  for(k in seq(K)){
    y_arr = y_arr + ws[k]*dnorm(x_arr, alpha_arr[k], beta_arr[k])
  }
  y_arr = y_arr + 1/L
  return(y_arr)
}

###### main simulation code ######
# mu_f = 350
# sigma_f = 30
# L = 1500
# 
# n_frag = 1000
# f_len_arr = round(rnorm(n_frag, mean=mu_f,sd = sigma_f ))
# 
# alpha_arr = sort(c(500, 800, 1000))
# n_ats = length(alpha_arr)
# beta_arr = sample(c(10,20,30), n_ats, replace = T)
# 
# unif_ws = 0.05
# ws = 1+runif(n_ats)
# ws = ws/sum(ws)
# ws = (1-unif_ws)*ws
#   
# boarder = round(n_frag*cumsum(ws))
# seg_st = c(1,boarder)
# seg_en = c(boarder-1,n_frag)
# 
# label_arr = rep(0,n_frag)
# st_arr = rep(0,n_frag)
# en_arr = rep(0,n_frag)
# for(i in seq(n_ats+1)){
#   
#   tmpinds = seg_st[i]:seg_en[i] 
#   tmpn = length(tmpinds)
#   
#   label_arr[tmpinds]=i
#   
#   if(i<=n_ats){ # ATS component
#     st_arr[tmpinds] = round(rnorm(tmpn, alpha_arr[i], beta_arr[i]))
#     en_arr[tmpinds] = st_arr[tmpinds] - f_len_arr[tmpinds]
#   }else{ # uniform component
#     st_arr[tmpinds] = sample(seq(L),tmpn,replace = T)
#     en_arr[tmpinds] = sample(seq(L),tmpn,replace = T)
#   }
# }
# 
# 
# # pair end
# res = atsmix(n_max_ats = 5, n_min_ats=1, st_arr = st_arr , en_arr=en_arr, L = L,
#              mu_f=350, # fragment length mean
#              sigma_f=30, # fragment length standard deviation
# 
#              # pa site information
#              min_ws = 0.01, # minimum weight of ATS site
# 
#              # inference with fixed parameters
#              fixed_inference_flag = FALSE,
# 
#              #single end mode
#              single_end_mode = FALSE,
#              debug=F)

# # single end
# res = atsmix(n_max_ats = 5, n_min_ats=1, st_arr = rep(NA,n_frag) , en_arr=en_arr, L = L,
#              mu_f=350, # fragment length mean
#              sigma_f=30, # fragment length standard deviation
#              
#              # pa site information
#              min_ws = 0.01, # minimum weight of ATS site
#              
#              # inference with fixed parameters
#              fixed_inference_flag = FALSE,
#              
#              #single end mode
#              single_end_mode = T, 
#              debug=F)
# 
# 
# cat('\n********************* Ground Truth **********************\n',sep='')
# cat("real ws: ",ws,"\n")
# cat("real alpha: ",alpha_arr,"\n")
# cat("real beta: ",beta_arr,"\n")
# 
# 
# plot(density(st_arr,bw='SJ'),col='red')
# lines(density(en_arr),col='blue')
# abline(v=res$alpha_arr,col='green')
# legend('topleft',legend=c('l','r','est_l'),lty=c(1,1,1),col=c('red','blue','green'))



# test_xu <- read.delim("~/Downloads/test_xu.txt")
# len_list = vector('list',length=nrow(test_xu))
# for(i in seq(nrow(test_xu))){
#   st_arr = as.numeric(strsplit(test_xu[i,'st_arr'],',')[[1]])
#   en_arr = as.numeric(strsplit(test_xu[i,'en_arr'],',')[[1]])
#   len_list[[i]] = st_arr-en_arr
# }
# len_arr = unlist(len_list)
# plot(density(len_arr),title('fragment size distribution'),xlim=c(0,500))


# for(i in c(6)){
#   st_arr = as.numeric(strsplit(test_xu[i,'st_arr'],',')[[1]])
#   en_arr = as.numeric(strsplit(test_xu[i,'en_arr'],',')[[1]])
#   L = max(st_arr+100)
  
  # res = atsmix(n_max_ats = 3, n_min_ats=2, st_arr = st_arr , en_arr=en_arr, L = L,
  #              mu_f=350, # fragment length mean
  #              sigma_f=30, # fragment length standard deviation
  #              
  #              # pa site information
  #              min_ws = 0.01, # minimum weight of ATS site
  #              
  #              # inference with fixed parameters
  #              fixed_inference_flag = FALSE,
  #              
  #              #single end mode
  #              single_end_mode = FALSE,
  #              debug=F)
  # 
  # 
  # 
  # plot(density(st_arr,bw='SJ'),col='red')
  # lines(density(en_arr),col='blue')
  # abline(v=res$alpha_arr,col='green')
  # legend('topleft',legend=c('l','r','est_l'),lty=c(1,1,1),col=c('red','blue','green'))
# }

# set.seed(5)
# # data = readRDS('1_975899-977241-.rds')
# data = readRDS('1_1212795-1214738-.rds')
# # data = readRDS('1_1307169-1308618-.rds')


# st_arr = abs(data$st_arr)
# en_arr = abs(data$en_arr)
# len_arr = st_arr-en_arr+1
# # plot(density(len_arr),title('fragment size distribution'))
# # plot(density(len_arr),title('fragment size distribution'),xlim=c(0,500))


# L = max(max(st_arr)+100, 2000)

# res = atsmix(n_max_ats = 3, n_min_ats=1, st_arr = st_arr , en_arr=en_arr, L = L,
#              # mu_f=300, # fragment length mean
#              mu_f = round(mean(len_arr)),
#              sigma_f=50, # fragment length standard deviation

#              # pa site information
#              min_ws = 0.01, # minimum weight of ATS site
#              max_unif_ws = 0.1,
#              max_beta = 50,

#              # inference with fixed parameters
#              fixed_inference_flag = FALSE,

#              #single end mode
#              single_end_mode = FALSE,
#              debug=F)


# min_val = min(min(st_arr), min(en_arr))
# max_val = max(max(st_arr), max(en_arr))
# plot(density(-st_arr),col='red',xlim=c(-max_val-100,0))
# lines(density(-en_arr),col='blue')
# # plot(density(-st_arr,bw='SJ'),col='red',xlim=c(-max_val-100,0))
# # lines(density(-en_arr,bw='SJ'),col='blue')
# x_arr = seq(-max_val-100,0)
# y_arr = est_density(-res$alpha_arr, res$beta_arr, res$ws, L, x_arr)
# lines(x_arr, y_arr, col='green')
# abline(v=-res$alpha_arr,col='green',lty=3)
# legend('topright',legend=c('l','r','est_l'),lty=c(1,1,1),col=c('red','blue','green'))

# # title( sub = paste0( 'ws=',paste(round(res$ws,digits = 2),collapse=" "),
# #                      ' mu_len=',round(mean(len_arr)), ' std_len=', round(sd(len_arr)))  )

# disp_ws = round(res$ws,digits = 2)
# K = length(disp_ws)-1
# disp_ws = c(disp_ws[K:1], disp_ws[K+1])
# title( sub = paste0( 'ws=',paste(disp_ws,collapse=" "),
#                      ' mu_len=',round(mean(len_arr)), ' std_len=', round(sd(len_arr)))  )





