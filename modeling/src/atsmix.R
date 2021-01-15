# ats mixture model inference
atsmix = function(n_max_ats, # maximum number of ATS sites
                n_min_ats=1, # minimum number of ATS sites
                
                # information for all DNA fragments, must be n_frag x 1 vector if specified
                st_arr, # l, start location of each DNA fragment, from 5' to 3'on the 5'-UTR
                en_arr, # r, end location of each DNA fragment

                L=NULL,  # length of UTR region
                
                # fragment size information
                mu_f=350, # fragment length mean
                sigma_f=50, # fragment length standard deviation
                
                # pa site information
                min_ws = 0.01, # minimum weight of ATS site
                max_beta = 50, # maximum std for ATS site
                
                # inference with fixed parameters
                fixed_inference_flag = FALSE,
                
                #single end mode
                single_end_mode = FALSE, 
                debug=T){

cal_z_k = function(alpha_arr, beta_arr, ws, k, log_zmat){
    K = length(ws) - 1  # last component is uniform component
    if(k<=K){
        log_zmat[,k] = log(ws[k]) + lik_lr_ab(st_arr, en_arr, alpha_arr[k], beta_arr[k], log=T)
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

mstep = function(alpha_arr, beta_arr, ws, Z, k){
    
    K = length(ws) - 1 # last component is uniform component
    new_ws = maximize_ws(Z)
    new_alpha_arr = alpha_arr
    new_alpha_arr[k] = sum(Z[,k]*st_arr)/sum(Z[,k])
    
    new_beta_arr = beta_arr;
    new_beta_arr[k] = sqrt( sum(Z[,k]*(st_arr-new_alpha_arr[k])^2)/sum(Z[,k]) )
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
em_algo = function(ws, para_mat, debug=F){
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
        temp = log_zmat
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

init_para = function(n_ats, debug=F){
    min_pos = min(st_arr)
    max_pos = max(st_arr)
    alpha_arr = min_pos + (max_pos-min_pos)*runif(n_ats)
    beta_arr = runif(n_ats,10,70)
    return(cbind(alpha_arr, beta_arr))
}

em_optim0 = function(n_ats, debug=F){
    n_trial=20
    
    lb_arr = rep(-Inf,n_trial)
    bic_arr = rep(Inf,n_trial)
    res_list = vector("list",n_trial)
    
    for(i in seq(n_trial)){                    
        res_list[[i]] = tryCatch({
            ws = runif(n_ats+1)+1  # last component is for uniform component
            ws[n_ats+1] = ws[n_ats+1]-1
            ws = ws/sum(ws)
            
            # initilize alpha_arr and beta_arr, considered pre_alpha_arr and pre_beta_arr
            para_mat = init_para(n_ats)
            
            em_algo(ws, para_mat, debug=debug)
        },error=function(e){
            print(e)
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
    res$beta_arr = round(res$beta_arr)
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

if (typeof(st_arr) == "list") {
    st_arr = unlist(st_arr)
}

if (typeof(en_arr) == "list") {
    en_arr = unlist(en_arr)
}

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
    res_list[[i]] = em_optim0(i, debug=debug)
    lb_arr[i] = tail(res_list[[i]]$lb_arr,1)
    bic_arr[i] = res_list[[i]]$bic
}

min_ind = which.min(bic_arr)
res = res_list[[min_ind]]

if(is.null(res$ws)){
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

return(res)
}