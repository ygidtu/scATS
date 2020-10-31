# apa mixture model inference
apamix = function(n_max_apa, # maximum number of APA sites
                  n_min_apa=1, # minimum number of APA sites
                  r1_utr_st_arr, # start location of each read on UTR part
                  r1_len_arr, # length of each read on UTR part
                  r2_len_arr, # length of each read on polyA part
                  polya_len_arr, # length of polyA
                  ws=NULL, # weights of each PA sites in the mixture
                  pre_theta_arr=NULL,   # pre-defined locations of PA sites
                  pre_theta_cnt_arr=NULL, # counts of reads corresponds to each predefined theta
                  L=NULL,  # length of UTR region
                  min_LA = 20, # minimum polyA length
                  max_LA=250, # maximum polyA length
                  mu_f=300, # fragment length mean
                  sigma_f=50, # fragment length standard deviation 
                  min_ws = 0.01, # minimum weight of PA site
                  min_gap_pa_sites = 50,  # minimum length between two PA sites
                  debug=T){
  
  estep = function(ws, x_arr, l_arr, r_arr, s_arr,  theta_arr){
    K = length(ws)
    N = length(x_arr)
    Z = matrix(0,nrow=N,ncol=K)
    for(k in seq(K)){
      Z[,k] = ws[k] * lik_f(x_arr, l_arr, r_arr, s_arr,  theta_arr[k])
    }
    return(Z/apply(Z,1,sum))
  }
  
  # calculate the value of kth PA site
  cal_theta = function(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr, k){
    K = length(ws)
    min_theta = ifelse(k==1, min_theta, theta_arr[k-1]+1)
    max_theta = ifelse(k==K, L, theta_arr[k+1]-1)
    all_theta_arr = seq(min_theta, max_theta, 10)
    
    #remove theta that are within 10bp of neighbouring PA sites
    tmp_gap=9
    for(tmp in c(theta_arr[k-1],theta_arr[k+1])){  # a=c(1,2) a[0]==numeric(0)  a[3]=NA
      if(is.na(tmp)) next()
      tmpinds = all_theta_arr>(tmp+tmp_gap) | all_theta_arr<(tmp-tmp_gap)
      all_theta_arr = all_theta_arr[tmpinds]
    }
    
    n_theta = length(all_theta_arr)
    if(n_theta==0){
      return(theta_arr[k])
    }
    
    lik_arr = rep(-Inf,n_theta)
    for(i in seq(n_theta)){
      tmp_theta_arr = theta_arr
      tmp_theta_arr[k] = all_theta_arr[i]
      # some reads may get NaN probability if this condition does not hold
      if(any(tmp_theta_arr>=min_max_theta)){
        lik_arr[i] = exp_log_lik(Z, ws, x_arr, l_arr, r_arr, s_arr,  tmp_theta_arr)
      }
    }
    max_ind = which.max(lik_arr)
    
    if(lik_arr[max_ind]==-Inf){
      return(theta_arr[k])
    }
    
    return(all_theta_arr[max_ind])
  }
  
  mstep = function(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr, k){
    N = length(x_arr)
    K = length(ws)
    
    new_ws = apply(Z,2,sum)/N
    new_theta_arr = theta_arr
    # when k=0, do not calculate theta_arr
    if(k>0){
      new_theta_arr[k] = cal_theta(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr, k)
    }
    return(list(ws=new_ws, theta_arr=new_theta_arr))
  }
  
  elbo = function(Z, ws, x_arr, l_arr, r_arr, s_arr, theta_arr){
    LZ = Z
    LZ[Z!=0] = log(Z[Z!=0])
    entropy = -1 * Z * LZ
    
    lb = exp_log_lik(Z, ws, x_arr, l_arr, r_arr, s_arr, theta_arr)+sum(entropy)
    if(is.na(lb)){
      stop("lower bounder is na.")
    }
    return (lb)
  }
  
  # calculate the expected log joint likelihood
  exp_log_lik = function(Z, ws, x_arr, l_arr, r_arr, s_arr, theta_arr){
    if(any(is.na(theta_arr))){
      print('some theta is NA')
    }
    
    K = length(ws)
    N = length(x_arr)
    ZZ = matrix(0,nrow=N,ncol=K)
    for(k in seq(K)){
      ZZ[,k] = log(ws[k]) + log(lik_f(x_arr, l_arr, r_arr, s_arr, theta_arr[k]))
    }
    ZZ[Z==0] = 0
    if(is.nan(sum(Z*ZZ))){
      stop('ZZ contains NaN.')
    }
    return(sum(Z*ZZ))
  }
  
  # calculate sum_s(f)
  lik_f = function(x_arr, l_arr, r_arr, s_arr, theta){
    ps = lik_s(s_dis_arr)   # s_dis_arr is in outer function
    ns = length(s_dis_arr)
    res = rep(0,length(x_arr))
    
    valid_inds = !is.na(s_arr)
    oth_inds = !valid_inds
    
    # calculate for read pairs with known s
    tmp1 = lik_l_xt(l_arr[valid_inds], x_arr[valid_inds], theta)
    tmp2 = lik_x_st(x_arr[valid_inds], s_arr[valid_inds], theta)
    res[valid_inds] = tmp1 * tmp2
    
    if( !any(oth_inds) ){
      return(res)
    }
    
    for(i in seq(ns)){
      s = s_dis_arr[i]  
      tmp1 = lik_r_s(r_arr[oth_inds], s)
      tmp2 = lik_x_st(x_arr[oth_inds], s, theta)
      tmp3 = lik_l_xt(l_arr[oth_inds], x_arr[oth_inds], theta)
      res[oth_inds] = res[oth_inds] + tmp1*tmp2*tmp3*ps[i]
      
      if(is.nan(sum(res))){
        stop('res contains NaN.')
      }
      
      if(any(is.infinite( res ))){
        stop('res contains Inf.')
      }
      
    }
    return(res)
  }
  
  # calculate p(l|x,theta)
  # l_arr is vector, x_arr is vector, theta is scalar
  lik_l_xt = function(l_arr, x_arr, theta){
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
    
    return( res )
  }
  
  # calculate p(x|sn,theta), x_arr is a vector, s and theta are scalar
  # s also be a vector when used as s_arr
  lik_x_st = function(x_arr, s, theta){
    res = dnorm(x_arr, theta+s+1-mu_f, sigma_f)
    if(any(is.infinite( res ))){
      stop('res contains Inf.')
    }
    return(res)
  }
  
  # calculate p(r|s), r is a vector, s is a scalar
  lik_r_s = function(r_arr,s){
    res = (r_arr<=s)/s
    if(any(is.infinite( res ))){
      stop('res contains Inf.')
    }
    return(res)
  }
  
  # calculate p(s)
  lik_s = function(s_dis_arr){
    # res = (s_dis_arr<=max_LA & s_dis_arr>=min_LA)/(max_LA-min_LA+1)
    n = length(s_dis_arr)
    res = rep(1/n, n)
    if(any(is.infinite( res ))){
      stop('res contains Inf.')
    }
    return( res )
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
  
  cal_bic = function(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr, pre_theta_arr=NULL){
    N = dim(Z)[1]
    K = length(ws)
    K1 = K - length(pre_theta_arr) # length(NULL)==0
    
    res = -2*exp_log_lik(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr) + (K+K1) *log(N)  # the smaller bic, the better model
    return(res)
  }
  
  # perform inference for K components
  em_algo = function(ws, x_arr, l_arr, r_arr, s_arr,  theta_arr, pre_theta_arr=NULL, pre_theta_cnt_arr=NULL,debug=F){
    #Z = estep(ws, x_arr, l_arr, r_arr, s_arr,  theta_arr)
    lb = -Inf
    lb_arr = rep(NA,nround)
    N = length(x_arr)
    
    # pre_theta_inds = which(theta_arr %in% pre_theta_arr) # index of predefined theta in theta_arr
    if(is.null(pre_theta_arr)){
      k_arr = gen_k_arr(length(theta_arr), nround)
      pre_cnt_arr = NULL
    }else{
      oth_theta_inds = which( !(theta_arr %in% pre_theta_arr) )  
      tmp_arr = gen_k_arr(length(oth_theta_inds), nround)
      k_arr = oth_theta_inds[tmp_arr]
      pre_cnt_arr = rep(0,length(theta_arr))
      pre_cnt_arr[theta_arr %in% pre_theta_arr] = pre_theta_cnt_arr
    }
    
    # in case all theta are given
    if(!is.null(pre_theta_arr) && length(oth_theta_inds)==0){
      Z = estep(ws, x_arr, l_arr, r_arr, s_arr,  theta_arr)
      res = mstep(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr, 0)
      pre_ws = (apply(Z,2,sum)+pre_cnt_arr)/(N+sum(pre_cnt_arr))  # ws plus pre_theta counts
      lb_new = elbo(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr)
      bic = cal_bic(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr, pre_theta_arr)
      if(debug){
        cat("estimated ws: ",res$ws,"\n")
        cat("estimated ws (plus pre_cnt): ",pre_ws,"\n")
        cat("estimated theta: ",theta_arr,"\n")
      }
      
      return(list(ws=ws, pre_ws=pre_ws, theta_arr=theta_arr, lb_arr=lb_new, bic=bic))
    }
    
    for(i in seq(nround)){
      if(debug){
        cat('iteration=',i,'  lb=',lb, "\n")
      }
      Z = estep(ws, x_arr, l_arr, r_arr, s_arr,  theta_arr)
      res = mstep(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr, k_arr[i])
      theta_arr = res$theta_arr
      ws = res$ws
      
      lb_new = elbo(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr)
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
        cat('Run all ',i,' iterations.',' lb=',lb, "\n")
      }
    }else{
      if(debug){
        cat('Converge in ',i,' iterations.',' lb=',lb, "\n")
      }
    }
    
    if(is.null(pre_theta_arr)){
      pre_ws = ws 
    }else{
      pre_ws = (apply(Z,2,sum)+pre_cnt_arr)/(N+sum(pre_cnt_arr)) 
    }
    
    bic = cal_bic(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr, pre_theta_arr)
    
    if(debug){
      cat("bic=",bic,"\n",sep = "")
      cat("estimated ws: ",ws,"\n")
      if(!is.null(pre_theta_arr)){
        cat("estimated ws (plus pre_cnt): ",pre_ws,"\n")
      }
      cat("estimated theta: ",theta_arr,"\n")
    }
    
    lb_arr = lb_arr[!is.na(lb_arr)]
    return(list(ws=ws, pre_ws=pre_ws, theta_arr=theta_arr, lb_arr=lb_arr, bic=bic))
  }
  
  # initialize theta
  # theta_arr is a vector of predifined theta that needs to be included
  init_theta = function(x_arr, l_arr, L, K, pre_theta_arr=NULL){
    min_theta = min(l_arr)
    min_max_theta = max(x_arr+l_arr-1)
    max_theta = L
    
    full_set = seq(min_theta,max_theta,10)
    
    if(!is.null(pre_theta_arr)){
      # check if input is valid
      stopifnot(!is.unsorted(pre_theta_arr))
      # pre_theta_arr = sort(pre_theta_arr)
      stopifnot(all(diff(pre_theta_arr)>=min_gap_pa_sites))
      tmp_k =length(pre_theta_arr)
      
      for(i in seq(tmp_k)){
        tmpinds = full_set<=(pre_theta_arr[i]-min_gap_pa_sites) | full_set>=(pre_theta_arr[i]+min_gap_pa_sites)
        full_set = full_set[tmpinds]
      }
    }else{
      tmp_k = 0
    }
    n_max_trial = 200
    for(i in seq(n_max_trial)){
      tmp_theta_arr = c( sample(full_set,K-tmp_k), pre_theta_arr)
      tmp_theta_arr = sort(tmp_theta_arr)
      flag1 = all( diff(tmp_theta_arr) >= min_gap_pa_sites)
      flag2 = tmp_theta_arr[K]>=min_max_theta
      if( flag1 & flag2 ){
        break
      }
    }
    if( i==n_max_trial ){
      stop(paste0('Failed to generate valid theta after ',n_max_trial,' trials.'))
    }
    
    return(tmp_theta_arr)
  }
  
  # merge thetas whose distance are smaller than the min_gap_pa_sites
  merge_theta = function(ws, theta_arr, pre_theta_arr=NULL){
    tmpind = which.min(diff(theta_arr))
    tmpind1 = tmpind+1
    stopifnot(min(diff(theta_arr))<min_gap_pa_sites)
    stopifnot(!all(theta_arr[c(tmpind,tmpind1)] %in% pre_theta_arr)) # two pre given pa sites cannot be the minimum
    
    if(theta_arr[tmpind] %in% pre_theta_arr){
      return( merge_comp(tmpind, tmpind1, ws, theta_arr) )
    }else if(theta_arr[tmpind1] %in% pre_theta_arr){
      return( merge_comp(tmpind1, tmpind, ws, theta_arr) )
    }else{
      return( merge_comp(tmpind1, tmpind, ws, theta_arr) )
    }
  }
  
  # merge component ind2 into component ind1
  merge_comp = function(ind1, ind2, ws, theta_arr){
    ws[ind1] = ws[ind2] + ws[ind1]
    ws = ws[-ind2]
    theta_arr = theta_arr[-ind2]
    return(list(ws=ws,theta_arr=theta_arr))
  }
  
  del_ws = function(ws, theta_arr, pre_theta_arr=NULL){
    if(length(theta_arr)==1){
      return( list(ws=ws, theta_arr=theta_arr) )
    }
    oth_theta_inds = which( !(theta_arr %in% pre_theta_arr) )
    tmpinds = which( ws[oth_theta_inds] < min_ws  )
    tmpinds = oth_theta_inds[tmpinds]
    
    tmpind = sample(tmpinds,1)
    ws = ws[-tmpind]
    ws = ws/sum(ws)
    theta_arr = theta_arr[-tmpind]
    
    return( list(ws=ws, theta_arr=theta_arr) )
  }
  
  em_optim0 = function(n_apa, L, x_arr, l_arr, r_arr, s_arr,  pre_theta_arr=NULL, pre_theta_cnt_arr=NULL, debug=F){
    n_trial=10
    
    lb_arr = rep(-Inf,n_trial)
    bic_arr = rep(Inf,n_trial)
    res_list = vector("list",n_trial)
    
    for(i in seq(n_trial)){
      if(debug){
        cat('\n-----------------K=',n_apa,' | ', 'i_trial=',i, ' | n_trial=',n_trial,' -------------\n',sep='')
      }
      
      res_list[[i]] = tryCatch({
        ws = runif(n_apa)+1
        ws = ws/sum(ws)
        
        theta_arr = init_theta(x_arr, l_arr, L, n_apa, pre_theta_arr)  # L and n_apa are from outer function
        # pre_theta_inds = which(theta_arr %in% pre_theta_arr) # index of predefined theta in theta_arr
        if(debug){
          cat("initial ws: ",ws,"\n")
          cat("initial theta: ",theta_arr,"\n")
        }
        
        em_algo(ws, x_arr, l_arr, r_arr, s_arr,  theta_arr, pre_theta_arr, pre_theta_cnt_arr, debug=debug)
      },error=function(e){
        list(ws=NULL,pre_ws=NULL,theta_arr=NULL,lb_arr=-Inf,bic=Inf)
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
      cat("theta: ",res$theta_arr,"\n")
      cat("lb=",tail(res$lb_arr,1),"\n")
      cat("bic=",tail(res$bic,1),"\n","\n")
    }else{
      cat("No results available for K=",n_apa,"\n\n")
    }
    
    return(res)
  }
  
  if(!is.null(pre_theta_arr)){
    stopifnot(n_max_apa>=length(pre_theta_arr))
    if(n_min_apa<=length(pre_theta_arr)){
      n_min_apa = length(pre_theta_arr)
    }
    
    tmp_idx = order(pre_theta_arr)
    pre_theta_arr = pre_theta_arr[tmp_idx]
    pre_theta_cnt_arr = pre_theta_cnt_arr[tmp_idx]
  }
  
  min_theta = min(r1_len_arr)
  min_max_theta = max(r1_utr_st_arr+r1_len_arr-1)
  
  n_frag = length(r1_len_arr)
  s_dis_arr = seq(min_LA,max_LA,10)  # 改成1
  stopifnot(n_frag==length(r1_utr_st_arr))
  stopifnot(n_frag==length(r2_len_arr))
  nround=200
  
  lb_arr = rep(-Inf,n_max_apa)
  bic_arr = rep(Inf,n_max_apa)
  ws_flag_arr = rep(FALSE, n_max_apa)
  res_list = vector("list",n_max_apa)
  
  for(i in seq(n_max_apa,n_min_apa)){
    cat('\n********************* K=',i, ' **********************\n',sep='')
    res_list[[i]] = em_optim0(i, L, r1_utr_st_arr, r1_len_arr, r2_len_arr, polya_len_arr, pre_theta_arr, pre_theta_cnt_arr,debug=debug)
    lb_arr[i] = tail(res_list[[i]]$lb_arr,1)
    bic_arr[i] = res_list[[i]]$bic
    ws_flag_arr[i] = all( res_list[[i]]$ws >= min_ws)
  }
  
  valid_inds = which(ws_flag_arr)
  
  min_ind = which.min(bic_arr[valid_inds])
  min_ind = valid_inds[min_ind]
  res = res_list[[min_ind]]
  
  if(debug){
    cat(paste0('ws=',sprintf("%.3f", res$ws),collapse=" "),"\n")
    idx = seq(3,nround)
    if(!is.na(res$lb_arr[idx[1]])){
      plot(idx,res$lb_arr[idx],type='o')
    }
  }
  
  cat('\n********************* Final Result **********************\n',sep='')
  cat("estimated ws: ",res$ws,"\n")
  cat("estimated theta: ",res$theta_arr,"\n")
  cat("lb=",tail(res$lb_arr,1),"\n")
  cat("bic=",res$bic,"\n")
  
  # calculate read assignment
  Z = estep(res$ws, r1_utr_st_arr, r1_len_arr, r2_len_arr, polya_len_arr, res$theta_arr)
  label = apply(Z, 1, which.max)
  res$label = label
  return(res)
}


set.seed(1)
new_s = rep(NA,length(s))
fix_s_inds = sample(length(s),round(0.5*length(s)))
new_s[fix_s_inds] = s[fix_s_inds]

# res = apamix(4, 1, x, l, r, L=1000, pre_theta_arr = pre_theta_arr, pre_theta_cnt_arr = pre_theta_cnt_arr,debug=T)
# res = apamix(4, 1, x, l, r, L=L, debug=F)
res = apamix(3, 2, x, l, r, new_s, L=1000,debug=F)

cat("     real ws: ",ws,"\n")
cat("     real theta: ",theta_arr,"\n")


# speical for real data
# n_frag = length(x)
# s = rep(100,n_frag) 
# LA = 250

coverage_cnt = rep(0,L+LA)
for(i in seq(n_frag)){
  coverage_cnt[x[i]:(x[i]+l[i]-1)] = coverage_cnt[x[i]:(x[i]+l[i]-1)] + 1
  coverage_cnt[(L+s[i]-r[i]+1):(L+s[i])] = coverage_cnt[(L+s[i]-r[i]+1):(L+s[i])] + 1
}

plot(seq(L+LA),coverage_cnt,type="s",
     # xlab = paste("UTR | polyA","ws=", paste(sprintf("%.2f",ws), collapse=" ")), 
     xlab = paste0("UTR | polyA"), 
     ylab = "coverage count",
     main = "red: PA sites, green: estimated PA sites, blue: UTR&polyA boarder")
abline(v=L, col="blue",lwd=3)
for(i in theta_arr){
  abline(v=i, col="red",lwd=2,lty=2)
}
for(i in res$theta_arr){
  abline(v=i, col="green",lwd=2,lty=2)
}