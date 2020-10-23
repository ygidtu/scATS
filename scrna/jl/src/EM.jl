

module EM
    using Distributions
    using Formatting
    using Memento
    using Random

    using Memento
    using Compat: @__MODULE__  # requires a minimum of Compat 0.26. Not required on Julia 0.7

    # Create our module level LOGGER (this will get precompiled)
    Memento.config!("debug"; fmt="[{date} - {level} | {name}]: {msg}")
    const LOGGER = getlogger(@__MODULE__)

    # Register the module level LOGGER at runtime so that folks can access the LOGGER via `get_LOGGER(MyModule)`
    # NOTE: If this line is not included then the precompiled `MyModule.LOGGER` won't be registered at runtime.
    function __init__()
        Memento.register(LOGGER)
    end
    
    Number = Union{Int64, Float64, Int}

    struct EMData
        ws::AbstractArray
        alpha_arr::Vector
        beta_arr::Vector
        lb_arr::Vector
        bic::Number
    end

    struct WinK
        alpha::Number
        beta::Number
        loglik::Number
    end

    export emptyData, lik_f0, em_algo

    function emptyData()::EMData
        return EMData(nothing, nothing, nothing, nothing, nothing)
    end

    function lik_f0(L::Number, max_LA::Number, do_log::Bool=true)::Float64
        pl_s = px = 1 / L
        pr_s = 1 / max_LA
        res = px * pl_s * pr_s

        if do_log
            return log(res)
        end
        return res
    end

    function lik_r_s(r_arr, s, do_log=false)
        res = [r_arr .<= s] ./ s
        if do_log
            return log.(res)
        else
            return res
        end
    end

    function lik_x_st(x_arr, s, theta, mu_f, sigma_f, do_log::Bool=false)
        if do_log
            return logpdf.(Normal(theta+s+1-mu_f, sigma_f), x_arr)
        else
            return pdf.(Normal(theta+s+1-mu_f, sigma_f), x_arr)
        end
    end
    
    function lik_l_xt(x_arr, l_arr, theta, log::Bool=false)
        utr_len_arr = theta .- x_arr .+ 1
        valid_inds = l_arr .<= utr_len_arr

        if any(isnan(utr_len_arr[valid_inds]))
            warn(LOGGER, "some length is 0.")
        end

        res[valid_inds] = 1 ./ utr_len_arr[valid_inds]
        if any(isinf(res))
            throw(string("res contains inf: ", res))
        end

        if log
            return log.(res)
        else
            return res
        end
    end

    function lik_lsr_t(
        x_arr, l_arr, r_arr, s_arr,
        theta, mu_f, sigma_f, pmf_s_dis_arr, 
        s_dis_arr, log::Bool = false)

        s_len = length(s_dis_arr)
        res = zeros(length(x_arr))

        oth_inds = isnan.(x_arr)
        valid_inds =  1 .- x_arr

        n_valid_inds = sum(valid_inds)
        n_other_inds = sum(oth_inds)

        if n_valid_inds > 0
            tmp1 = lik_x_st(x_arr[valid_inds], s_arr[valid_inds], theta, mu_f, sigma_f, log)
            tmp2 = lik_l_xt(x_arr[valid_inds], s_arr[valid_inds], theta, log)

            if log
                res[valid_inds] = tmp1 .+ tmp2
            else
                res[valid_inds] = tmp1 .* tmp2
            end
        end

        if !any(oth_inds)
            return res
        end

        tmp_sum = zeros(n_other_inds)
        for i = 1:length(s_len)
            tmp1 = lik_r_s(r_arr[oth_inds], s_dis_arr[i])
            tmp2 = lik_x_st(x_arr[oth_inds], s_dis_arr[i], theta, mu_f, sigma_f, log)
            tmp3 = lik_l_xt(x_arr[oth_inds], l_arr[oth_inds], theta, log)
            tmp_sum .+= tmp1 .* tmp2 .* tmp3 .* pmf_s_dis_arr[i]
        end
        
        if log
            res[oth_inds] = log.(tmp_sum)
        else
            res[oth_inds] = tmp_sum
        end
        return res
    end

    function lik_lsr_t0(x_arr, l_arr, r_arr, s_arr, theta, pmf_s_dis_arr, s_dis_arr, log)
        s_len = length(s_dis_arr)

        try
            res = zeros(len(x_arr))
            oth_inds = isnan.(s_arr)
            valid_inds = 1 .- oth_inds
        catch e
            x_arr, l_arr, r_arr, res = [x_arr], [l_arr], [r_arr], zeros(1)
            oth_inds = [isnan.(s_arr)]
            valid_inds = [1 .- oth_inds]
        end

        n_valid_inds = sum(valid_inds)
        n_other_inds = sum(oth_inds)

        if n_valid_inds > 0
            res[valid_inds] = lik_l_xt(x_arr[valid_inds], l_arr[valid_inds], theta, log)
        end
        if ! any(oth_inds)
            return res
        end

        tmp_sum = zeros(n_other_inds)
        for i = 1:length(s_len)
            tmp1 = lik_r_s(r_arr[oth_inds], s_dis_arr[i])
            tmp2 = 1
            tmp3 = lik_l_xt(x_ar[oth_inds], l_arr[oth_inds], theta, log)

            tmp_sum .+= tmp1 .* tmp2 .* tmp3 .* pmf_s_dis_arr[i]
        
        end

        if log
            res[oth_inds] = log(tmp_sum)
        else
            res[oth_inds] = tmp_sum
        end
        return res
    end

    function maximize_ws(Z)
        return sum(Z, dims = 1) ./ size(Z)[1]
    end

    function maximize_win_k(
        theta_inds::Int, beta_ind::Int, zk::Vector,
        all_theta::Vector, init_beta_arr::Vector, all_win_log_lik_mat_list::Vector
    )::WinK
        n_input_theta = length(theta_inds)
        valid_inds = findall(x -> x > 0, zk)
        
        beta = init_beta_arr[beta_ind]  # init_beta_arr defined in outer function
        all_win_log_lik_mat = all_win_log_lik_mat_list[beta_ind]
        
        tmp_lik_arr = zeros(n_input_theta)
        for i in 1:n_input_theta
          tmpvec = all_win_log_lik_mat[valid_inds, theta_inds[i]]
          tmp_lik_arr[i]= sum(zk[valid_inds] .* tmpvec)
        end
  
        max_ind = findfirst(x -> x == maximum(tmp_lik_arr), tmp_lik_arr)
        max_loglik = tmp_lik_arr[max_ind]
        alpha = all_theta[theta_inds[max_ind]]
        return(WinK(alpha, beta, max_loglik))
    end

    function exp_log_lik(log_zmat, Z)::Number
        ZZ = Z .* log_zmat
        ZZ[findall(isequal(0), Z)] .= 0
        return sum(ZZ)
    end

    function elbo(log_zmat, Z)
        LZ = deepcopy(Z)
        LZ[findall(!isequal(0), Z)] .= log.(Z[findall(!isequal(0), Z)])
        entropy = -1 .* Z .* LZ
        lb = exp_log_lik(log_zmat, Z) .+ sum(entropy)
        if isnan(lb)
            throw(string("lower bounder is na: ", lb))
        end
        return lb
    end

    function cal_z_k(
        all_win_log_lik_mat_list::Vector,
        alpha_arr::Vector, beta_arr::Vector, ws, k::Int,
        log_zmat, all_theta::Vector, init_beta_arr::Vector,
        unif_log_lik::Number)

        K = length(ws) - 1
        tmp_len = length(all_win_log_lik_mat_list)
        # print(beta_arr)
        if k <= K
            tmp_ind = [x for x in 1:tmp_len][init_beta_arr .== beta_arr[k]]
            all_win_log_lik_mat = all_win_log_lik_mat_list[tmp_ind[1]]
            log_zmat[:, k] = log(ws[k]) .+ collect(Iterators.flatten(all_win_log_lik_mat[:, all_theta .== alpha_arr[k]]))
        else
            log_zmat[:, K+1] .= log(ws[K+1]) + unif_log_lik
        end
        return log_zmat
    end

    function norm_z(Z::AbstractArray)
        Z = Z .- maximum(Z, dims=2)
        Z = [exp.(x) for x in Z]
        return Z ./ sum(Z, dims=2)
    end

    function cal_bic(log_zmat, Z)
        N, K = size(Z)
        K = K - 1
        return -2 * exp_log_lik(log_zmat, Z) + (3 * K + 1) * log(N)
    end

    function fixed_inference(
        alpha_arr::Vector, beta_arr::Vector, n_frag::Number, nround::Int, 
        all_theta::Vector, all_win_log_lik_mat_list::Vector, init_beta_arr::Vector,
        unif_log_lik::Number, verbose::Bool=false)

        lb, lb_arr = -Inf, Array{Number}(undef, nround, 1)

        K = length(alpha_arr)

        ws = [rand(Uniform()) + 1 for _ in 1:K+1]
        ws = ws .+ 1
        ws = ws ./ sum(ws)

        k_arr = gen_k_arr(K, nround)
        log_zmat = zeros(n_frag, K + 1)

        for k = 1:K
            log_zmat = cal_z_k(
                all_win_log_lik_mat_list,
                alpha_arr, beta_arr, ws, k,
                log_zmat, all_theta, init_beta_arr, unif_log_lik
            )
        end

        Z = nothing
        i = 1
        for i = 1:nround
            if verbose
                info(LOGGER, Formatting.format(
                    FormatExpr("iteration={}, lb={}"), i, lb
                ))
            end

            log_zmat = cal_z_k(
                all_win_log_lik_mat_list,
                alpha_arr, beta_arr, ws,
                k_arr[i], log_zmat, all_theta,
                init_beta_arr, unif_log_lik
            )
            
            Z = norm_z(log_zmat)
            ws = maximize_ws(Z)

            lb_new = elbo(log_zmat, Z)
            lb_arr[i] = lb_new
            if lb_new == -Inf
                lb = -Inf
                break
            end
            if abs(lb_new - lb)  < abs(1e6 * lb)
                break
            else
                lb = lb_new
            end
        end
        if verbose
            info(LOGGER, Formatting.format(
                FormatExpr("Run all {} iterations. lb={}"), 
                i, lb
            ))
        end
        bic = NaN
        if !isnothing(Z)
            bic = cal_bic(log_zmat, Z)
        end
        
        # label = argmax(log_zmat, dims=2)
        if verbose
            info(LOGGER, Formatting.format(
                FormatExpr("bic={}; estimated ws: {}; estimated alpha:  {}; estimated beta: {}"), 
                bic, ws, alpha_arr, beta_arr
            ))
        end

        lb_arr = lb_arr[findall(x -> !isnan(x), lb_arr)]
        return EMData(ws, alpha_arr, beta_arr, lb_arr, bic)
    end

    function mstep(
        alpha_arr::Vector, beta_arr::Vector, ws::AbstractArray, 
        Z::AbstractArray, k::Int, all_theta::Vector, min_theta::Number,
        L::Int, init_beta_arr::Vector, all_win_log_lik_mat_list::Vector;
        pre_init_alpha_arr=nothing
    )::EMData
        K = length(ws) - 1
        new_ws = maximize_ws(Z)

        if !isnothing(pre_init_alpha_arr) && sum([x in pre_init_alpha_arr for x in alpha_arr[k]]) > 0
            valid_inds = findall(x -> x == all_theta, alpha_arr[k])[1]
        else
            if k == 1 || k > K
                tmp_min_theta = min_theta
            else
                tmp_min_theta = alpha_arr[k-1] + beta_arr[k-1]
            end

            if k >= K
                tmp_max_theta = L
            else
                tmp_max_theta = alpha_arr[k + 1]
            end

            valid_inds = findall(x -> x == 1, [tmp_min_theta <= i <= tmp_max_theta for i = all_theta])[1]
        end

        n_winsize = length(init_beta_arr)
        res_list = Vector{WinK}()
        loglik_arr = [-Inf for _ in 1:n_winsize]

        for i in 1:n_winsize
            push!(
                res_list, 
                maximize_win_k(
                    valid_inds, i, Z[:, k], 
                    all_theta, init_beta_arr, all_win_log_lik_mat_list
                )
            )
            loglik_arr[i] = res_list[i].loglik
        end
        res = res_list[findfirst(isequal(maximum(loglik_arr)), loglik_arr)]
    
        alpha_arr[k] = res.alpha
        beta_arr[k] = res.beta

        return EMData(new_ws, alpha_arr, beta_arr, Vector(), NaN)
    end

    function gen_k_arr(K::Int64, n::Int64)::Vector{Int64}

        if K <= 1
            return [K for _ in 1:n]
        else
            res = Vector()
            count_index = 0
            pool = [x for x in 1:K]
            last = nothing
            while count_index < n
                count_index += 1
                pool = shuffle(pool)
    
                if pool[1] == last
                    swap_with = rand(1:length(pool))
                    pool[1], pool[swap_with] = pool[swap_with], pool[1]
                end
    
                append!(res, pool)
                last = pool[length(pool)]
            end
    
            return res
        end
    end

    function em_algo(
        ws::AbstractArray, theta_win_mat, n_frag::Int, pre_alpha_arr::Vector,
        pre_beta_arr::Vector, all_theta::Vector, init_beta_arr::Vector,
        all_win_log_lik_mat_list::Vector, pre_init_alpha_arr::Vector,
        min_pa_gap::Int, unif_log_lik::Number, min_theta::Number,
        L::Int; nround::Int=50, verbose::Bool=false)::EMData

        # n_frag=self.n_frag
        # pre_alpha_arr=self.pre_alpha_arr
        # pre_beta_arr=self.pre_beta_arr
        # all_theta=self.all_theta
        # init_beta_arr=self.init_beta_arr
        # all_win_log_lik_mat_list=self.all_win_log_lik_mat_list
        # pre_init_alpha_arr=self.pre_init_alpha_arr
        # min_pa_gap=self.min_pa_gap
        # unif_log_lik=self.unif_log_lik
        # min_theta=self.min_theta
        # L=self.utr_len

        lb = -Inf
        lb_arr = [NaN for _ = 1:nround]
        K = length(ws) - 1

        # assume this is a NamedTuple
        alpha_arr = theta_win_mat.res_alpha_arr
        beta_arr = theta_win_mat.res_beta_arr

        if length(pre_alpha_arr) > 0
            pre_alpha_arr = [i in pre_init_alpha_arr for i in alpha_arr]
        end
        if length(pre_alpha_arr) > 0 && length(pre_beta_arr) > 0
            tmp_pre_init_alpha_arr = pre_init_alpha_arr[1 .- isnan.(pre_beta_arr)]
            tmp_arr = gen_k_arr(length(alpha_arr) - length(tmp_pre_init_alpha_arr), nround)
            opt_inds = [1:length(alpha_arr)][1 .- [i in tmp_pre_init_alpha_arr for i in alpha_arr]]
            k_arr = tmp_arr

            for i = 1:length(opt_inds)
                k_arr[tmp_arr == i] = opt_inds[i]
            end
        elseif length(pre_alpha_arr) > 0 && len(pre_beta_arr) < 1
            k_arr = gen_k_arr(length(alpha_arr), nround)
        elseif length(pre_alpha_arr) < 1
            pre_init_alpha_arr = nothing
            k_arr = gen_k_arr(length(alpha_arr), nround)
        end
     
        log_zmat = zeros(n_frag, K+1)
       
        for k = 1:(K+1)
            log_zmat = cal_z_k(
                all_win_log_lik_mat_list,
                alpha_arr, beta_arr,
                ws, k, log_zmat,
                all_theta, init_beta_arr, unif_log_lik
            )
        end

        Z = nothing
        i = 1
        for i = 1:nround
            log_zmat = cal_z_k(
                all_win_log_lik_mat_list,
                alpha_arr, beta_arr,
                ws, k_arr[i],
                log_zmat, all_theta,
                init_beta_arr, unif_log_lik
            )
   
            Z = norm_z(log_zmat)

            res = mstep(
                alpha_arr, beta_arr,
                ws, Z, k_arr[i],
                all_theta, min_theta, L,
                init_beta_arr,
                all_win_log_lik_mat_list,
                pre_init_alpha_arr=pre_init_alpha_arr
            )
    
            alpha_arr = res.alpha_arr
            beta_arr = res.beta_arr
            ws = res.ws

            lb_new = elbo(log_zmat, Z)
            lb_arr[i] = lb_new
            if lb_new == -Inf
                lb = -Inf
                break
            end
            
            if abs(lb_new - lb) < abs(1e-6 * lb)
                break
            else
                lb = lb_new
            end
        end

        if verbose
            info(LOGGER, Formatting.format(
                FormatExpr("Run all {} interactions. lb={}"), i, lb
            ))
        end

        bic = NaN
        if !isnothing(Z) 
            bic = cal_bic(log_zmat, Z)
        end
        
        if verbose
            info(LOGGER, Formatting.format(
                FormatExpr("bic={}; estimated ws: {}; estimated alpha:  {}; estimated beta: {}"), 
                bic, ws, alpha_arr, beta_arr
            ))
        end

        if length(pre_alpha_arr) > 0
            oth_inds = findall(isequal(0), pre_flag_arr)

            for i in oth_inds
                if any(abs(alpha_arr[i] .- pre_alpha_arr) .< min_pa_gap)
                    warn(LOGGER, Formatting.format(
                        FormatExpr("alpha: {} is within {} distance"),
                        alpha_arr[i], min_pa_gap, 

                    ))
                    return emptyData()
                end
            end
        end

        lb_arr = lb_arr[findall(x -> !isnan(x), lb_arr)]
        return EMData(ws, alpha_arr, beta_arr, lb_arr, bic)
    end
end

