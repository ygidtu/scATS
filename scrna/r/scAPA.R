# 66e01bd
# apa mixture model inference
# apa mixture model inference
apamix = function(n_max_apa, # maximum number of APA sites
                  n_min_apa=1, # minimum number of APA sites
                  
                  # information for all reads, must be n_frag x 1 vector if specified
                  r1_utr_st_arr, # start location of each read on UTR part of this fragment, x
                  r1_len_arr, # length of each read on UTR part, l
                  r2_len_arr, # length of each read on polyA part, r
                  polya_len_arr=NULL, # length of polyA of this fragment, s
                  pa_site_arr=NULL, # pa site locations of reads, theta
                                    # reads cover junction sites have value, use NA for other reads
                  
                  L=NULL,  # length of UTR region
                  
                  # polya length information
                  empirical_LA_flag = TRUE, # use empirical polya distribution or not
                  LA_dis_arr = c(10,30,50,70,90,110,130), # empirical polya lengths
                  pmf_LA_dis_arr = c(309912, 4107929, 802856, 518229, 188316, 263208, 101), # empirical pmf of corresponding polya lengths
                  min_LA=20, # minimum polyA length
                  max_LA=150, # maximum polyA length
                  
                  # fragment size information
                  mu_f=270, # fragment length mean
                  sigma_f=30, # fragment length standard deviation
                  
                  # pa site information
                  min_ws = 0.01, # minimum weight of PA site
                  min_pa_gap = 100, # minimum distance between PA sites
                  max_beta = 70, # maximum standard deviation of PA sites
                  
                  # pre-specified PA sites
                  pre_alpha_arr = NULL,  # pre-specified alpha values, e.g. c(500, 1000)
                  pre_beta_arr = NULL,   # pre-specified beta values. e.g. c(40,NA). If all beta values are unknown, set it to NULL.
                  
                  # inference with fixed parameters
                  fixed_inference_flag = FALSE,
                  pre_init_theta_win_mat = NULL, # fixed initialization parameters, e.g. cbind(c(500,1150,1300), c(70,10,10)) represent pre_alpha_arr and pre_beta_arr
                  debug=T){
  
  cal_z_k = function(alpha_arr, beta_arr, ws, k, log_zmat){
    K = length(ws) - 1  # last component is uniform component
    if(k<=K){
      all_win_log_lik_mat = all_win_log_lik_mat_list[[which(init_beta_arr==beta_arr[k])]]
      log_zmat[,k] = log(ws[k]) + all_win_log_lik_mat[, which(all_theta==alpha_arr[k])]
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
    return( ws )
  }
  
  # maximize alpha_k given beta
  # choose alpha and beta of window k that maximize the likelihood
  # theta_inds: index of theta in all_theta the locate between component k-1 and k+1
  # n_step_win: number of steps (theta_step) this window contains
  # zk: the kth column of Z, assignment probability of each read to kth component
  maximize_win_k = function(theta_inds, beta_ind, zk){
    n_input_theta = length(theta_inds)
    valid_inds = which(zk>0)
    
    beta = init_beta_arr[beta_ind]  # init_beta_arr defined in outer function
    all_win_log_lik_mat = all_win_log_lik_mat_list[[beta_ind]]
    
    tmp_lik_arr = rep(0, n_input_theta)
    for(i in seq(n_input_theta)){
      tmpvec = all_win_log_lik_mat[valid_inds,theta_inds[i]]
      tmp_lik_arr[i]= sum(zk[valid_inds] * tmpvec)
    }
    max_ind = which.max(tmp_lik_arr)
    max_loglik = tmp_lik_arr[max_ind]
    alpha = all_theta[theta_inds[max_ind]]
    return(list(alpha=alpha, beta=beta, loglik=max_loglik))
  }
  
  mstep = function(alpha_arr, beta_arr, ws, Z, k, pre_init_alpha_arr=NULL){
    K = length(ws) - 1 # last component is uniform component
    new_ws = maximize_ws(Z)
    
    if(alpha_arr[k] %in% pre_init_alpha_arr){
      valid_inds = which(all_theta==alpha_arr[k]) # no need to update alpha if given
    }else{
      tmp_min_theta = ifelse(k==1, min_theta, alpha_arr[k-1]+beta_arr[k-1])   # min_theta defined in outer function
      tmp_max_theta = ifelse(k==K, L, alpha_arr[k+1])
      valid_inds = which( all_theta>=tmp_min_theta & all_theta<=tmp_max_theta ) # all_theta defined in outer function
    }
    
    n_winsize = length(init_beta_arr)
    res_list = vector("list", length=n_winsize)
    loglik_arr = rep(-Inf, n_winsize)

    for(i in c(seq(n_winsize))){
      res_list[[i]] = maximize_win_k(valid_inds, i, Z[,k])
      loglik_arr[i] = res_list[[i]]$loglik
    }
    max_ind = which.max(loglik_arr)
    res = res_list[[max_ind]]
    
    alpha_arr[k] = res$alpha
    beta_arr[k] = res$beta
    
    return(list(alpha_arr=alpha_arr, beta_arr=beta_arr, ws=new_ws))
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
  
  # calculate likelihood for the uniform component, it is a constant number
  lik_f0 = function(log=F){
    # ps = pmf_s_dis_arr # pmf_s_dis_arr is in outer function
    px = 1/L
    pl_s = 1/L
    pr_s = 1/max_LA   # max_LA defined in outer function
    res = px*pl_s*pr_s
    
    if(log){
      return(log(res))
    }else{
      return(res)
    }
  }
  
  # calculate p(x,l,r|theta)  or  log ( p(x,l,r|theta) )
  lik_lsr_t = function(x_arr, l_arr, r_arr, s_arr, theta, log=F){
    ps = pmf_s_dis_arr # pmf_s_dis_arr is in outer function
    ns = length(s_dis_arr)
    res = rep(0,length(x_arr))
    
    # s_arr are length of polya part for each pair-end reads, 
    valid_inds = !is.na(s_arr)  # pair end reads which touch PA sites, i.e. known polya length for these reads
    oth_inds = !valid_inds
    
    n_valid_inds = sum(valid_inds)
    n_other_inds = sum(oth_inds)
    
    if(n_valid_inds>0){
      tmp1 = lik_l_xt(l_arr[valid_inds], x_arr[valid_inds], theta, log=log)
      tmp2 = lik_x_st(x_arr[valid_inds], s_arr[valid_inds], theta, log=log)
      if(log){
        res[valid_inds] = tmp1 + tmp2
      }else{
        res[valid_inds] = tmp1 * tmp2
      }
    }
    
    if( !any(oth_inds) ){
      return(res)
    }
    
    tmp_sum = rep(0,n_other_inds)
    for(i in seq(ns)){
      s = s_dis_arr[i]
      tmp1 = lik_r_s(r_arr[oth_inds], s)
      tmp2 = lik_x_st(x_arr[oth_inds], s, theta)
      tmp3 = lik_l_xt(l_arr[oth_inds], x_arr[oth_inds], theta)
      tmp_sum = tmp_sum + tmp1*tmp2*tmp3*ps[i]
    }
    if(log){
      res[oth_inds] = log(tmp_sum)
    }else{
      res[oth_inds] = tmp_sum
    }
    
    return(res)
  }
  
  # calculate p(x,l,r|theta)  or  log ( p(x,l,r|theta) ), where theta is known
  lik_lsr_t0 = function(x_arr, l_arr, r_arr, s_arr, theta, log=F){
    ps = pmf_s_dis_arr # pmf_s_dis_arr is in outer function
    ns = length(s_dis_arr)
    res = rep(0,length(x_arr))
    
    # s_arr are length of polya part for each pair-end reads, 
    valid_inds = !is.na(s_arr)  # pair end reads which touch PA sites, i.e. known polya length for these reads
    oth_inds = !valid_inds
    
    n_valid_inds = sum(valid_inds)
    n_other_inds = sum(oth_inds)
    
    if(n_valid_inds>0){
      tmp1 = lik_l_xt(l_arr[valid_inds], x_arr[valid_inds], theta, log=log)
      
      # tmp2 = lik_x_st(x_arr[valid_inds], s_arr[valid_inds], theta, log=log)
      if(log){
        tmp2 = 0
        res[valid_inds] = tmp1 + tmp2
      }else{
        tmp2 = 1
        res[valid_inds] = tmp1 * tmp2
      }
    }
    
    if( !any(oth_inds) ){
      return(res)
    }
    
    tmp_sum = rep(0,n_other_inds)
    for(i in seq(ns)){
      s = s_dis_arr[i]
      tmp1 = lik_r_s(r_arr[oth_inds], s)
      tmp2 = 1
      # tmp2 = lik_x_st(x_arr[oth_inds], s, theta)
      tmp3 = lik_l_xt(l_arr[oth_inds], x_arr[oth_inds], theta)
      tmp_sum = tmp_sum + tmp1*tmp2*tmp3*ps[i]
    }
    if(log){
      res[oth_inds] = log(tmp_sum)
    }else{
      res[oth_inds] = tmp_sum
    }
    
    return(res)
  }
  
  # calculate p(l|x,theta)
  # l_arr is vector, x_arr is vector, theta is scalar
  lik_l_xt = function(l_arr, x_arr, theta, log=F){
    utr_len_arr = theta - x_arr + 1
    valid_inds = (l_arr <= utr_len_arr)
    res = valid_inds+0
    if(any(is.na(utr_len_arr[valid_inds]))){
      print("some length is 0.")
    }
    res[valid_inds] = 1/utr_len_arr[valid_inds]
    
    if(any(is.infinite( res ))){
      stop('res contains Inf.')
    }
    if(log){
      return( log(res) )
    }else{
      return( res )
    }
  }
  
  # calculate p(x|sn,theta), x_arr is a vector, s and theta are scalar
  # s also be a vector when used as s_arr
  lik_x_st = function(x_arr, s, theta, log=F){
    res = dnorm(x_arr, theta+s+1-mu_f, sigma_f, log = log)
    return(res)
  }
  
  # calculate p(r|s), r is a vector, s is a scalar
  lik_r_s = function(r_arr,s, log=F){
    res = (r_arr<=s)/s
    if(log){
      return( log(res) )
    }else{
      return(res)
    }
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
  
  # replace query with closest values in ref_arr
  replace_with_closest = function(ref_arr, query_arr){
    n = length(query_arr)
    res = rep(0,n)
    for(i in seq(n)){
      tmpind = which.min( abs(ref_arr-query_arr[i]) )
      res[i] = ref_arr[tmpind]
    }
    return(res)
  }
  
  # perform inference given alpha_arr and beta_arr
  fixed_inference = function(alpha_arr, beta_arr, debug=F){
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
      if(debug){
        cat('iteration=',i,'  lb=',lb, "\n")
      }
      
      # estep
      log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k_arr[i], log_zmat)
      Z = norm_z(log_zmat)
      
      # res = mstep(alpha_arr, beta_arr, ws, Z, k_arr[i]) # no need to 
      
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
    if(i==nround){
      if(debug){
        cat('Run all ',i,' iterations.\n','lb=',lb, "\n")
      }
    }else{
      if(debug){
        cat('Converge in ',i,' iterations.\n','lb=',lb, "\n")
      }
    }
    
    bic = cal_bic(log_zmat, Z)
    label = apply(log_zmat, 1, which.max)
    
    if(debug){
      cat("bic=",bic,"\n",sep = "")
      cat("estimated ws: ",ws,"\n")
      cat("estimated alpha: ",alpha_arr,"\n")
      cat("estimated beta: ",beta_arr,"\n")
    }
    
    lb_arr = lb_arr[!is.na(lb_arr)]
    if(debug){
      nd = length(lb_arr)
      if(nd>=3){
        plot(seq(nd-2),lb_arr[3:nd],type='o')
      }
    }
    return(list(ws=ws, alpha_arr=alpha_arr, beta_arr=beta_arr, lb_arr=lb_arr, bic=bic, label=label))
  }
  
  # perform inference for K components
  # theta_arr => theta_win_mat
  em_algo = function(ws, theta_win_mat, debug=F){
    lb = -Inf
    lb_arr = rep(NA,nround)
    N = n_frag
    K = length(ws)-1
    
    alpha_arr = theta_win_mat[,1]
    beta_arr = theta_win_mat[,2]
    
    if(!is.null(pre_alpha_arr)){
      pre_flag_arr = alpha_arr %in% pre_init_alpha_arr
    }
    
    if(!is.null(pre_alpha_arr) && !is.null(pre_beta_arr)){
      tmp_pre_init_alpha_arr = pre_init_alpha_arr[!is.na(pre_beta_arr)]  # pre_init_alpha_arr defined in outer function
      k_arr = gen_k_arr(length(alpha_arr)-length(tmp_pre_init_alpha_arr), nround)
      opt_inds = which(!alpha_arr %in% tmp_pre_init_alpha_arr)
      tmp_arr = k_arr
      for(i in seq(length(opt_inds))){
        tmp_arr[k_arr==i] = opt_inds[i]
      }
      k_arr = tmp_arr
    }else if(!is.null(pre_alpha_arr) && is.null(pre_beta_arr)){
      k_arr = gen_k_arr(length(alpha_arr), nround)
    }else if(is.null(pre_alpha_arr)){
      pre_init_alpha_arr = NULL
      k_arr = gen_k_arr(length(alpha_arr), nround)
    }
    
    log_zmat = matrix( 0, nrow=N, ncol=K+1 )
    for(k in seq(K+1)){
      log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k, log_zmat)
    }
    
    for(i in seq(nround)){
      if(debug){
        cat('iteration=',i,'  lb=',lb, "\n")
      }
      
      # estep
      log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k_arr[i], log_zmat)
      Z = norm_z(log_zmat)

      res = mstep(alpha_arr, beta_arr, ws, Z, k_arr[i], pre_init_alpha_arr)

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
    if(i==nround){
      if(debug){
        cat('Run all ',i,' iterations.\n','lb=',lb, "\n")
      }
    }else{
      if(debug){
        cat('Converge in ',i,' iterations.\n','lb=',lb, "\n")
      }
    }
    
    bic = cal_bic(log_zmat, Z)
    
    if(debug){
      cat("bic=",bic,"\n",sep = "")
      cat("estimated ws: ",ws,"\n")
      cat("estimated alpha: ",alpha_arr,"\n")
      cat("estimated beta: ",beta_arr,"\n")
    }
    
    if(!is.null(pre_alpha_arr)){
      other_inds = which(!pre_flag_arr)  # note estimated alpha may equal pre_alpha, there could be two alpha with the same value, thus need to use pre_flag_arr computed in the beginning
      for(i in other_inds){
        if( any( abs(res$alpha_arr[i]-pre_alpha_arr) < min_pa_gap ) ){
          warning(paste0('alpha[',i,']=',res$alpha_arr[i],' is within min_pa_gap=',min_pa_gap, 'bp of pre_alpha_arr=',
                         paste0(pre_alpha_arr,collapse = " ")))
          return(list(ws=NULL,alpha_arr=NULL,beta_arr=NULL,lb_arr=-Inf,bic=Inf))
        }
      }
    }
    
    lb_arr = lb_arr[!is.na(lb_arr)]
    if(debug){
      nd = length(lb_arr)
      if(nd>=3){
        plot(seq(nd-2),lb_arr[3:nd],type='o')
      }
    }
    return(list(ws=ws, alpha_arr=alpha_arr, beta_arr=beta_arr, lb_arr=lb_arr, bic=bic))
  }
  
  # calculate coverage 
  split_data = function(x_arr, l_arr, L){
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
    
    n_frag = length(x_arr)
    coverage_cnt = rep(0,L)
    for(i in seq(n_frag)){
      tmpinds = x_arr[i]-1+seq(l_arr[i])
      coverage_cnt[tmpinds] = coverage_cnt[tmpinds] + 1
    }
    # kernel smoothing of coverage
    ks_res = ksmooth(seq(L), coverage_cnt, kernel='normal', bandwidth = 5*theta_step, x.points = seq(L)) # theta_step defined in outer function
    # change points
    sign_arr = sign(diff(c(-1,ks_res$y)))
    sign_arr = handle_zero(sign_arr)
    
    chng = cumsum( rle(sign_arr)$lengths )
    mode_arr = ks_res$x[ chng[seq(1,length(chng),2)] ]
    n_mode = length(mode_arr)
    boarder_arr = ks_res$x[ chng[seq(2,length(chng),2)] ][1:(n_mode-1)]
    st_arr = c(1,boarder_arr+1)
    en_arr = c(boarder_arr, L)
    ws_arr = rep(1,n_mode)
    for(i in seq(n_mode)){
      ws_arr[i] = sum(coverage_cnt[st_arr[i]:en_arr[i]])
    }
    ws_arr = ws_arr/sum(ws_arr)
    
    return(list(n_win=n_mode, st_arr=st_arr, en_arr=en_arr, ws_arr=ws_arr, mode_arr=mode_arr))
  }
  
  init_theta_win = function(K){
    sample_theta = function(k_arr, mode='full'){
      res_arr = rep(0,length(k_arr))
      for(k in seq(length(k_arr))){
        iw = k_arr[k]
        if(mode=='full'){
          tmp_theta_arr = all_theta[all_theta>=data_wins$st_arr[iw] & all_theta<data_wins$en_arr[iw]]
        }else if(mode=='left'){
          tmp_theta_arr = all_theta[all_theta>=data_wins$st_arr[iw] & all_theta<data_wins$mode_arr[iw]]
        }else if(mode=='right'){
          tmp_theta_arr = all_theta[all_theta>=data_wins$mode_arr[iw] & all_theta<data_wins$en_arr[iw]]
        }else{
          stop(paste0('Unknown mode: ',mode, 'Must be [left|right|full]'))
        }
        if(length(tmp_theta_arr)>0){
          res_arr[k] = sample(tmp_theta_arr,1)
        }else{
          res_arr[k] = -1
        }
      }
      return(sort(res_arr))
    }
    # special treatment for the second sampling
    sample_theta1 = function(k_arr){
      tmpinds = which( k_arr==(n_data_win+1) )
      if(length(tmpinds)>0){
        th1_arr = sample_theta(n_data_win, mode='full')
      }else{
        th1_arr = integer(0)
      }
      tmpinds = which( k_arr!=(n_data_win+1) )
      if(length(tmpinds)>0){
        th2_arr = sample_theta(k_arr[tmpinds], mode='left')
      }else{
        th2_arr = integer(0)
      }
      return( c(th1_arr, th2_arr) )
    }
    
    try_sample = function(){
      if(n_data_win>=K){
        k_arr = sort( sample(n_data_win, K, prob=data_wins$ws_arr) )
        return( sample_theta(k_arr, mode='right') )
      }else if(n_data_win>=K-K1){
        th_arr1 = sample_theta(seq(n_data_win), mode='right')
        if(K1==1){
          k_arr = k_big_inds+1
        }else{
          k_arr = sample(k_big_inds+1, K-n_data_win, prob=data_wins$ws_arr[k_big_inds])
        }
        th_arr2 = sample_theta1(k_arr)
        return( sort( c(th_arr1, th_arr2) ) )
      }else{
        th_arr1 = sample_theta(seq(n_data_win), mode='right')
        th_arr2 = sample_theta1(k_big_inds+1)
        k_arr = sample(n_data_win, K-n_data_win-K1, prob=data_wins$ws_arr, replace=T)
        th_arr3 = sample_theta(k_arr, mode='full')
        return( sort( c(th_arr1, th_arr2, th_arr3) ) )
      }
    }
    
    # process pre_alpha_arr and pre_beta_arr
    pre_para_flag = !is.null(pre_alpha_arr) # if true, pre_init_alpha_arr and pre_init_beta_arr are provided by outer function
    if(pre_para_flag){
      pre_init_alpha_arr = replace_with_closest(all_theta, pre_alpha_arr)
      pre_init_beta_arr = sample(init_beta_arr[seq(length(init_beta_arr)/2+1)], length(pre_alpha_arr), replace=T)
      if(!is.null(pre_beta_arr)){
        tmp_inds = which(!is.na(pre_beta_arr))
        pre_init_beta_arr[tmp_inds] = replace_with_closest(init_beta_arr, pre_beta_arr[tmp_inds])
      }
    }
    
    if(pre_para_flag && length(pre_alpha_arr)==K){
      return(cbind(pre_init_alpha_arr, pre_init_beta_arr))
    }
    
    n_data_win = data_wins$n_win  # data_win_list defined in outer function
    k_big_inds = which(data_wins$ws_arr>0.1) # big peaks
    K1 = length(k_big_inds)
    
    n_max_trial = 200
    for(i in seq(n_max_trial)){
      tmp_alpha_arr = try_sample()
      flag1 = all( diff(tmp_alpha_arr) >= min_pa_gap )
      flag2 = all( tmp_alpha_arr>0 )
      if(flag1 & flag2){
        break
      }
    }
    
    if( i==n_max_trial ){
      stop(paste0('Failed to generate valid theta after ',n_max_trial,' trials.'))
    }
    
    get_data_win_lab = function(alpha_arr){
      n_data_win = data_wins$n_win
      lab_vec = alpha_arr
      left = 0
      right = data_wins$mode_arr[1]
      for(i in seq(n_data_win)){
        lab_vec[alpha_arr>=left & alpha_arr<right] = i
        left = data_wins$mode_arr[i]
        right = data_wins$mode_arr[i+1]
      }
      lab_vec[alpha_arr>=left] = n_data_win+1
      return(lab_vec)
    }
    # replace part of sampled tmp_alpha_arr with pre_alpha_arr
    if(pre_para_flag){
      # get the component id of data_win for sampled alpha and beta
      sam_lab = get_data_win_lab(tmp_alpha_arr)
      pre_lab = get_data_win_lab(pre_alpha_arr)
      n_to_rm = 0
      for(i in seq(length(pre_lab))){
        smp_inds = which(sam_lab==pre_lab[i] & sam_lab>0)
        if(length(smp_inds)==0){
          n_to_rm = n_to_rm + 1
        }else if(length(smp_inds)==1){
          tmpind = smp_inds
          sam_lab[tmpind] = 0
        }else{
          tmpind = which.min( abs(tmp_alpha_arr[smp_inds]-pre_alpha_arr[i]) )
          sam_lab[smp_inds[tmpind]] = 0
        }
      }
      if(n_to_rm>0){
        smp_inds = which(sam_lab>0) # should have at least two items
        stopifnot(length(smp_inds)>1)
        tmpinds = sample(smp_inds,n_to_rm)
        sam_lab[tmpinds] = 0
      }
      valid_inds = which(sam_lab>0)
      res_alpha_arr = sort(c(tmp_alpha_arr[valid_inds], pre_init_alpha_arr))
      # find index of pre_alpha_arr
      pre_inds = match(pre_init_alpha_arr, res_alpha_arr)
      res_beta_arr = sample(init_beta_arr[seq(length(init_beta_arr)/2+1)],K, replace=T)
      res_beta_arr[pre_inds] = pre_init_beta_arr
    }else{
      res_alpha_arr = tmp_alpha_arr
      res_beta_arr = sample(init_beta_arr[seq(length(init_beta_arr)/2+1)],K, replace=T)
    }
    
    return(cbind(res_alpha_arr, res_beta_arr)) # return a matrix of two columns, K x 2 matrix
  }
  
  kernel_smooth = function(all_theta, all_theta_lik_mat, beta){
    n_step = beta/theta_step   # theta_step defined in outer function
    n_bd_step = ceiling(n_step * 3)      # 3 std away
    weights = dnorm(seq(-n_bd_step,n_bd_step), 0, n_step) 
    weights = weights/sum(weights)
    wlen = length(weights)
    
    n_all_theta = length(all_theta)
    n_frag = dim(all_theta_lik_mat)[1]
    stopifnot(n_all_theta==dim(all_theta_lik_mat)[2])
    
    all_win_log_lik_mat = matrix(0, nrow=n_frag, ncol=n_all_theta)
    
    for(i in seq(n_bd_step)){
      tmpinds = seq(i-n_bd_step,i+n_bd_step)
      tmpinds = tmpinds[tmpinds>0]
      tmpws = weights[tmpinds]
      tmpws = tmpws/sum(tmpws)
      all_win_log_lik_mat[,i] = log( apply( tmpws * t(all_theta_lik_mat[,tmpinds]) , 2, sum) )
    }
    
    for(i in seq(n_bd_step+1, n_all_theta-n_bd_step)){
      tmpinds = seq(i-n_bd_step,i+n_bd_step)
      tmpws = weights
      all_win_log_lik_mat[,i] = log( apply( tmpws * t(all_theta_lik_mat[,tmpinds]) , 2, sum) )
    }
    
    for(i in seq(n_all_theta-n_bd_step+1,n_all_theta)){
      tmpinds = seq(i-n_bd_step,i+n_bd_step)
      
      tmpws = weights[tmpinds<=n_all_theta]
      tmpws = tmpws/sum(tmpws)
      
      tmpinds = tmpinds[tmpinds<=n_all_theta]
      
      all_win_log_lik_mat[,i] = log( apply( tmpws * t(all_theta_lik_mat[,tmpinds]) , 2, sum) )
    }
    
    return(all_win_log_lik_mat)
  } 
  
  em_optim0 = function(n_apa, debug=F){
    n_trial=5
    
    lb_arr = rep(-Inf,n_trial)
    bic_arr = rep(Inf,n_trial)
    res_list = vector("list",n_trial)
    
    for(i in seq(n_trial)){
      if(debug){
        cat('\n-----------------K=',n_apa,' | ', 'i_trial=',i, ' | n_trial=',n_trial,' -------------\n',sep='')
      }
      
      res_list[[i]] = tryCatch({
        ws = runif(n_apa+1)+1  # last component is for uniform component
        ws[n_apa+1] = ws[n_apa+1]-1
        ws = ws/sum(ws)
        
        # initilize alpha_arr and beta_arr, considered pre_alpha_arr and pre_beta_arr
        theta_win_mat = init_theta_win(n_apa)
        
        if(debug){
          cat("initial ws: ",ws,"\n")
          cat("initial alpha: ",theta_win_mat[,1],"\n")
          cat("initial beta: ",theta_win_mat[,2],"\n")
        }

        em_algo(ws, theta_win_mat, debug=debug)
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
    
    cat('\n-----------------K=',n_apa,' final result -------------\n',sep='')
    if(!is.null(res$ws)){
      cat("ws: ",round(res$ws,digits=3),"\n")
      cat("alpha: ",res$alpha_arr,"\n")
      cat("beta: ",res$beta_arr,"\n")
      cat("lb=",tail(res$lb_arr,1),"\n")
      cat("bic=",tail(res$bic,1),"\n","\n")
    }else{
      cat("No results available for K=",n_apa,"\n\n")
    }
    
    return(res)
  }
  
  # remove components with weight less than min_ws
  rm_component = function(res, min_ws){
    K = length(res$alpha_arr)
    
    
    if(!is.null(pre_alpha_arr)){
      pre_flag = res$alpha_arr %in% pre_init_alpha_arr  # pre_init_alpha_arr defined in outer function
      rm_inds = which( (res$ws[1:K]<min_ws) & (!pre_flag) )
    }else{
      rm_inds = which(res$ws[1:K]<min_ws)
    }
    
    if(length(rm_inds)==0){
      return(res)
    }else{
      cat('Remove components ', rm_inds, ' with weight less than min_ws=',min_ws,"\n")
      cat(paste("ws=",round(res$ws[rm_inds],digits=3)), " " )
      cat(paste("alpha=",res$alpha_arr[rm_inds]), " ")
      cat(paste("beta=",res$beta_arr[rm_inds]), "\n")
      
      alpha_arr = res$alpha_arr[-rm_inds]
      beta_arr = res$beta_arr[-rm_inds]
      res = fixed_inference(alpha_arr,beta_arr)
      return(res)
    }
  }
  
  min_theta = min(r1_len_arr)
  if(!is.null(pre_alpha_arr)){
    min_theta = min(min_theta, min(pre_alpha_arr))
    stopifnot(max(pre_alpha_arr)<=L)
  }

  n_frag = length(r1_len_arr)
  if(empirical_LA_flag){
    s_dis_arr = LA_dis_arr #c(10,30,50,70,90,110,130)
    pmf_s_dis_arr = pmf_LA_dis_arr #c(309912, 4107929, 802856, 518229, 188316, 263208, 101)
  }else{
    s_dis_arr = seq(min_LA,max_LA,10)
    pmf_s_dis_arr = rep(1/length(s_dis_arr),length(s_dis_arr))
  }
  pmf_s_dis_arr = pmf_s_dis_arr/sum(pmf_s_dis_arr)
  stopifnot(n_frag==length(r1_utr_st_arr))
  stopifnot(n_frag==length(r2_len_arr))
  stopifnot(L>=min_theta)
  nround=50
  
  if(is.null(pa_site_arr)){
    pa_site_arr = rep(NA,n_frag)
  }
  
  if(is.null(polya_len_arr)){
    polya_len_arr = rep(NA,n_frag)
  }
  
  n_pre_pa_sites = 0
  if(!is.null(pre_alpha_arr)){
    if(is.unsorted(pre_alpha_arr)){
      stop(paste0('pre_alpha_arr=', paste0(pre_alpha_arr, collapse = " "), ' is not sorted.'))
    }
    if(max(pre_alpha_arr)>L){
      stop(paste0('All pre_alpha_arr=', paste0(pre_alpha_arr, collapse = " "), ' should be less than L=',L))
    }
    n_pre_pa_sites = length(pre_alpha_arr)
  }
  if(!is.null(pre_beta_arr)){
    if(n_pre_pa_sites!=length(pre_beta_arr)){
      stop( paste0('pre_alpha_arr=', paste0(pre_alpha_arr, collapse = " "), 
                   '  pre_beta_arr=', paste0(pre_beta_arr, collapse = " "), 
                   '  lengths are not the same') )
    }
    if(all(is.na(pre_beta_arr))){
      stop( paste0('All pre_alpha_arr=', paste0(pre_beta_arr, collapse = " "), ' are NA, set it to NULL instead.'))
    }
    if(max(pre_beta_arr, na.rm = T)>max_beta){
      stop(paste0('All pre_beta_arr=', paste0(pre_beta_arr, collapse = " "), ' should be less than max_beta=',max_beta))
    }
  }
  if(n_pre_pa_sites>0){
    if(n_pre_pa_sites>n_max_apa){
      warning(paste0("n_max_apa=",n_max_apa," changed to n_pre_pa_sites=",n_pre_pa_sites))
      n_max_apa = n_pre_pa_sites
    }
    if(n_pre_pa_sites>n_min_apa){
      warning(paste0("n_min_apa=",n_min_apa," changed to n_pre_pa_sites=",n_pre_pa_sites))
      n_min_apa = n_pre_pa_sites
    }
  }
  
  lb_arr = rep(-Inf,n_max_apa)
  bic_arr = rep(Inf,n_max_apa)
  ws_flag_arr = rep(FALSE, n_max_apa)
  res_list = vector("list",n_max_apa)
  
  theta_step = 9
  all_theta = seq(min_theta,L,theta_step)
  n_all_theta = length(all_theta)
  
  # calculate the likelihood for each theta
  all_theta_lik_mat = matrix(0, nrow=n_frag, ncol=n_all_theta)
  pa_inds = which(!is.na(pa_site_arr))
  oth_inds = which(is.na(pa_site_arr))
  for(i in seq(n_all_theta)){
    # all_theta_lik_mat[,i] = lik_lsr_t(r1_utr_st_arr, r1_len_arr, r2_len_arr, polya_len_arr, all_theta[i])
    all_theta_lik_mat[oth_inds,i] = lik_lsr_t(r1_utr_st_arr[oth_inds], r1_len_arr[oth_inds], r2_len_arr[oth_inds], polya_len_arr[oth_inds], all_theta[i])
  }
  for(i in pa_inds){
    tmpind = which.min(abs(all_theta-pa_site_arr[i]))
    # all_theta_lik_mat[i, tmpind] = lik_lsr_t(r1_utr_st_arr[i], r1_len_arr[i], r2_len_arr[i], polya_len_arr[i], all_theta[tmpind])
    all_theta_lik_mat[i, tmpind] = lik_lsr_t0(r1_utr_st_arr[i], r1_len_arr[i], r2_len_arr[i], polya_len_arr[i], all_theta[tmpind])
  }
  
  stopifnot(max_beta>=theta_step)
  init_beta_arr=seq(theta_step, ceiling(max_beta/theta_step)*theta_step,theta_step)
  # init_beta_arr = c(10,20,30,40,50,60,70)
  all_win_log_lik_mat_list = vector("list",length(init_beta_arr))
  for(i in seq(length(init_beta_arr))){
    all_win_log_lik_mat_list[[i]] = kernel_smooth(all_theta, all_theta_lik_mat, init_beta_arr[i] )
  }
  
  # process pre_alpha_arr and pre_beta_arr
  if(!is.null(pre_alpha_arr)){
    pre_init_alpha_arr = replace_with_closest(all_theta, pre_alpha_arr)
  }else{
    pre_init_alpha_arr = NULL
  }
  
  unif_log_lik = lik_f0(log = T)
  
  # provide results pre-specified parameters
  if(fixed_inference_flag && !is.null(pre_init_theta_win_mat)){
    alpha_arr = replace_with_closest(all_theta, pre_init_theta_win_mat[,1])  # derive approximate result
    beta_arr = replace_with_closest(init_beta_arr, pre_init_theta_win_mat[,2])
    res = fixed_inference(alpha_arr, beta_arr, debug = T)
    res$alpha_arr = pre_init_theta_win_mat[,1]
    res$beta_arr = pre_init_theta_win_mat[,2]
    return(res)
  }
  
  data_wins = split_data(r1_utr_st_arr, r1_len_arr, L)
  if(debug){
    cat('Initialization windows\n')
    tmp = cbind(data_wins$st_arr,data_wins$mode_arr, data_wins$en_arr, data_wins$ws_arr)
    colnames(tmp)=c('start','mode','end','weight')
    print(tmp)
  }

  for(i in seq(n_max_apa,n_min_apa)){
    cat('\n********************* K=',i, ' **********************\n',sep='')
    res_list[[i]] = em_optim0(i, debug=debug)
    lb_arr[i] = tail(res_list[[i]]$lb_arr,1)
    bic_arr[i] = res_list[[i]]$bic
  }

  min_ind = which.min(bic_arr)
  res = res_list[[min_ind]]
  
  if(is.null(res$ws)){
    cat('\n Inference failed. No results available. \n')
    return(res)
  }
  
  if(debug){
    cat(paste0('ws=',sprintf("%.3f", res$ws),collapse=" "),"\n")
    nd = length(res$lb_arr)
    if(nd>=3){
      plot(seq(nd-2),res$lb_arr[3:nd],type='o')
    }
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
  
  # replace pre_init_alpha_arr with original value, keep pre_init_beta_arr
  if(!is.null(pre_alpha_arr)){
    pre_inds = which(res$alpha_arr %in% pre_init_alpha_arr)  # potential problem, inferred alpha may be the same as pre_init_alpha_arr, leads to dupliates
    for(i in pre_inds){
      tmpind = which(res$alpha_arr[i]==pre_init_alpha_arr)
      res$alpha_arr[i] = pre_alpha_arr[tmpind]
    }
    
    # if(!is.null(pre_beta_arr)){
    #   n_pre = length(pre_inds)
    #   for(i in seq(n_pre)){
    #     if(!is.na(pre_beta_arr[i])){
    #       res$beta_arr[pre_inds[i]] = pre_beta_arr[i]
    #     }
    #   }
    #   }
  }
  
  cat('\n********************* Final Result **********************\n',sep='')
  cat("estimated ws: ",res$ws,"\n")
  cat("estimated alpha: ",res$alpha_arr,"\n")
  cat("estimated beta: ",res$beta_arr,"\n")
  cat("lb=",tail(res$lb_arr,1),"\n")
  cat("bic=",res$bic,"\n")
  
  return(res)
}

# run_expriment = function(df, filename){
#   offset = 150
#   df$V1 = df$V1 + offset
#   df$V7 = df$V7 + offset
  
#   L <- df[which.max(df$V1), 'V1'] + df[which.max(df$V2), 'V2'] + 50
  
#   res = apamix(
#     n_max_apa = 4,
#     n_min_apa = 2,
#     df$V1,
#     df$V2,
#     df$V3,
#     df$V6,
#     df$V7,
#     L = L,
#     mu_f=300,
#     sigma_f = 30,
#     pre_alpha_arr = c(20, 660),
#     # pre_beta_arr = c(30, NA),
#     # pre_init_theta_win_mat = cbind(c(500,1150,1300), c(70,10,10)),
#     debug=F
#   )
  
#   # res$alpha = res$alpha - offset
#   # df$V1 = df$V1 - offset
  
#   x = df$V1
#   l = df$V2
#   r = df$V3
#   s = df$V6
#   L = L
#   LA = 130
#   n_frag = dim(df)[1]
  
#   s[is.na(s)]=LA
#   coverage_cnt = rep(0,L+LA)
#   for(i in seq(n_frag)){
#     coverage_cnt[x[i]:(x[i]+l[i]-1)] = coverage_cnt[x[i]:(x[i]+l[i]-1)] + 1
#     coverage_cnt[(L+s[i]-r[i]):(L+s[i])] = coverage_cnt[(L+s[i]-r[i]):(L+s[i])] + 1
#   }
  
#   plot(seq(L+LA),coverage_cnt,type="s",
#        # xlab = paste("UTR | polyA","ws=", paste(sprintf("%.2f",ws), collapse=" ")), 
#        xlab = paste0("UTR | polyA"), 
#        ylab = "coverage count",
#        main = paste0(filename, " ", paste( round(res$ws,digits = 2), sep = " ", collapse = " ") ) ,
#        xaxt='n')
#   abline(v=L, col="green",lwd=3)
  
#   ymax = par("usr")[4]
  
#   for(i in seq(length(res$alpha_arr)) ){
#     st = res$alpha_arr[i] - res$beta_arr[i]
#     en = res$alpha_arr[i] + res$beta_arr[i]
#     x = c(st,st,en,en)
#     y = c(0,ymax,ymax,0)
#     lines(x,y,type='S',col='blue',lwd=2,lty=2)
#     abline(v=res$alpha_arr[i], col="yellow",lwd=2,lty=1)
#   }
  
#   axis(side = 1, at = c(seq(0,L+LA,100)), las=2 )
  
#   jrt = as.data.frame(table(df$V7[!is.na(df$V7)]))
#   loc = as.numeric(levels(jrt[,1]))
#   pmf = jrt[,2]/max(jrt[,2])
#   lines(loc, ymax*pmf, type='h', col='red', lwd=1)
  
#   # res$alpha_arr = res$alpha_arr - offset
#   return(res)
# }

# set.seed(1)

# ####### test data

# # pmf_polya_len = read.table(file="d_polya_len.txt", header = T)[,2]
# # pmf_polya_len = pmf_polya_len/sum(pmf_polya_len)

# # file = "/Users/lcheng/Documents/projects/scrna/zhouran/ENSMUSG00000000827.txt"
# # file = "/Users/lcheng/Documents/projects/scrna/zhouran/10_128799101_128800824.txt"
# # file = "/Users/lcheng/Documents/projects/scrna/zhouran/ENSMUSG00000000827_V7.txt"

# # file = "/Users/lcheng/Documents/projects/scrna/zhouran/problem_genes/1_15843083_15844052_Terf1.txt"

# input_dir = "/Users/lcheng/Documents/projects/scrna/zhouran/problem_genes/"
# file_list = list.files(input_dir)

# output_dir = "/Users/lcheng/Documents/projects/scrna/zhouran/problem_genes_figures/"

# # for(i in seq(1,length(file_list))){
# for(i in c(2)){
  
#   filename = file_list[i]
#   filestem = substr(filename,1,nchar(filename)-4)
#   file = paste0(input_dir,filename)
  
#   df <- read.table(file, row.names = NULL)
  
#   # df$V7 <- df$V1 + df$V2 -1
#   # df$V7[df$V5 == "no"] <- NA
  
#   cat( 'Processing ', i, 'th file: ', filename, '\n')
  
#   outfile = paste0(output_dir, filestem, '.pdf')
  
#   res = run_expriment(df, filestem)

#   # tryCatch({
#   #     pdf(file=outfile, width = 12, height = 6)
#   #     res = run_expriment(df, filestem)
#   #     dev.off()},
#   #          error = function(e){cat( 'error in ', i, 'th file: ', filename, '\n')})
# }





