# this script simulate isofroms and reads
source('sim_class.R')  # import all necessary functions
set.seed(5)

disp_info = function(i_frag){
  x = i_frag
  y = frag_label[x]
  cat(paste0("i_frag=",x,"\ti_isoform=",y,"\n"))
  cat(paste0("R1ST=",R1ST[x,y],"\tR1EN=", R1EN[x,y],"\n"))
  cat(paste0("R2ST=",R2ST[x,y],"\tR2EN=", R2EN[x,y],"\n"))
  cat(paste0("st_frag_on_gene=",st_frag_on_gene[x],"\ten_frag_on_gene=",en_frag_on_gene[x],"\n"))
  cat(paste0("local isoform i=",y,"\n"))
  cat(paste0('-----------------------\n'))
  print(local_iso_arr[[y]])
  cat(paste0("isoform i=",y," on gene\n"))
  cat(paste0('-----------------------\n'))
  print(iso_arr[[y]])
  cat(paste0("R1 and R2 of isoform i=",y," on gene\n"))
  cat(paste0('-----------------------\n'))
  map_winset(fromWinSet = local_iso_arr[[y]], 
             toWinSet = iso_arr[[y]], 
             subWinSet = winset_factory(win_list_factory(start_pos_arr = c(R1ST[x],R2ST[x]), end_pos_arr = c(R1EN[x], R2EN[x])) )
  )
}


construct_iso_arr = function(iso_on_gene){ # iso_on_gene is a list of K elements, each contains the start and end of exons, sorted according to start
  K = length(iso_on_gene)
  iso_arr = vector('list',K)
  local_iso_arr = vector('list',K)
  for(i in seq(K)){
    iso_arr[[i]] = winset_factory(win_list_factory(iso_on_gene[[i]]$start, iso_on_gene[[i]]$end))
    local_iso_arr[[i]] = get_local_winset(iso_arr[[i]])
  }
  return(list(local_iso_arr=local_iso_arr, iso_arr=iso_arr))
  
}



cal_relative_pos_on_isoform = function(n_frag, # 1x1, number of all fragment
                                       n_isoform, # 1x1, number of isoforms (K)
                                       frag_label, # nx1, label (which isoform) of each pair-end read
                                       R1ST, # n x K, start position of R1 on isoform k
                                       R1EN, # n x K, end position of R1 on isoform k
                                       R2ST, # n x K, start position of R2 on isoform k
                                       R2EN, # n x K, end position of R2 on isoform k
                                       local_iso_arr, 
                                       iso_arr){
  
  for(n in seq(n_frag)){
    if(n%%1000 == 0){
      cat(paste0('n=', n, " ", round(n*100/n_frag, 2), "%", "\n"))
    }
    
    # real_iso_ind = frag_label[n]
    real_iso_ind = 3
    n = 2449
    
    tmp1_st = R1ST[n, real_iso_ind]
    tmp1_en = R1EN[n, real_iso_ind]
    tmp2_st = R2ST[n, real_iso_ind]
    tmp2_en = R2EN[n, real_iso_ind]
    
    # 把R1 split后所从属的俩外显子的位置拿出来，R2所在的外显子的位置拿出来
    tmp1_winset = intersect(local_iso_arr[[real_iso_ind]],window_factory(tmp1_st,tmp1_en))
    tmp2_winset = intersect(local_iso_arr[[real_iso_ind]],window_factory(tmp2_st,tmp2_en))


    # fromWinSet = local_iso_arr[[real_iso_ind]]
    # toWinSet = iso_arr[[real_iso_ind]]
    # subWinSet = tmp1_winset

    winset = fromWinSet
    ref_ind = 915

    tmp1_read_gene_winset = map_winset(local_iso_arr[[real_iso_ind]], iso_arr[[real_iso_ind]], tmp1_winset) # read on gene
    tmp2_read_gene_winset = map_winset(local_iso_arr[[real_iso_ind]], iso_arr[[real_iso_ind]], tmp2_winset) # read on gene
    
    if(length(tmp1_winset@windows)>1){
      junc_read1_label[n] = T
    }
    if(length(tmp2_winset@windows)>1){
      junc_read2_label[n] = T
    }
    
    
    for(k in seq(n_isoform)){
      if(k==real_iso_ind){
        next
      }
      
      if(tmp1_read_gene_winset<=iso_arr[[k]] && tmp2_read_gene_winset<=iso_arr[[k]]){
        tmp1_read_local_k_winset = map_winset(iso_arr[[k]], local_iso_arr[[k]], tmp1_read_gene_winset)
        tmp2_read_local_k_winset = map_winset(iso_arr[[k]], local_iso_arr[[k]], tmp2_read_gene_winset)
        if(is_consecutive(tmp1_read_local_k_winset) && is_consecutive(tmp2_read_local_k_winset)){
          tmpst1 = get_start_pos(tmp1_read_local_k_winset)
          tmpen1 = get_end_pos(tmp1_read_local_k_winset)
          tmpst2 = get_start_pos(tmp2_read_local_k_winset)
          tmpen2 = get_end_pos(tmp2_read_local_k_winset)
          R1ST[n,k] = tmpst1[1]
          R1EN[n,k] = tail(tmpen1,1)
          
          R2ST[n,k] = tmpst2[1]
          R2EN[n,k] = tail(tmpen2,1)
          
          stopifnot(S[n,k]>0)
        }
      }
    }
  }
  return(list(rel_start_mat=R1ST, rel_end_mat=R2EN))
}


# sample code
load('data1.RData')
gtf_iso_arr = vector('list',n_isoform)
gtf_iso_arr[[1]] = list(start=c(1548, 2573, 4036), end=c(2032, 3063, 4539 ))
gtf_iso_arr[[2]] = list(start=c(2573, 4163), end=c(3063, 4652))
gtf_iso_arr[[3]] = list(start=c(1548, 2706, 4163 ), end=c(2032, 3198, 4652))

res = construct_iso_arr(gtf_iso_arr)
iso_arr = res$iso_arr
local_iso_arr = res$local_iso_arr

res_pos = cal_relative_pos_on_isoform(n_frag, # 1x1, number of all fragment
                                       n_isoform, # 1x1, number of isoforms (K)
                                       frag_label, # nx1, label (which isoform) of each pair-end read
                                       R1ST, # n x K, start position of R1 on isoform k
                                       R1EN, # n x K, end position of R1 on isoform k
                                       R2ST, # n x K, start position of R2 on isoform k
                                       R2EN, # n x K, end position of R2 on isoform k
                                       local_iso_arr, 
                                       iso_arr)

print(res_pos$rel_start_mat[3002,])
print(res_pos$rel_end_mat[3002,])
