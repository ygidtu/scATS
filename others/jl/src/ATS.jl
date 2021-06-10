module ATSMIX
    # using BSON
    using Distributions
    using Formatting
    using KernelDensity
    using Memento
    using OrderedCollections
    using Parameters
    using Random
    using RCall
    using Statistics
    using StatsBase

    include(joinpath(@__DIR__, "genomic.jl"))

    using Compat: @__MODULE__  # requires a minimum of Compat 0.26. Not required on Julia 0.7

    # Create our module level LOGGER (this will get precompiled)
    const LOGGER = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")

    # const LOGGER = getlogger(@__MODULE__)

    # Register the module level LOGGER at runtime so that folks can access the LOGGER via `get_LOGGER(MyModule)`
    # NOTE: If this line is not included then the precompiled `MyModule.LOGGER` won't be registered at runtime.
    function __init__()
        Memento.register(LOGGER)
    end

    Number = Union{Int64, Float64, Int}

    export emptyData, fit, fit_by_R, setLevel, EMHeader

    function setLevel(level::String)
        setlevel!(LOGGER, level)
    end

    @with_kw mutable struct Param
        st_arr::Vector
        en_arr::Vector
        mu_f::Number
        sigma_f::Number
        unif_log_lik::Number
        L::Int
        max_beta::Number
        seed::Int
        n_win::Number
        mode_arr::Vector
        step_size::Number
        fixed_inference_flag::Bool
        predef_beta_arr::Vector
        single_end_mode::Bool
        n_max_ats::Int = 5
        n_min_ats::Int = 1
        min_ws::Number = 0.01
        max_unif_ws::Number = 0.1
        cage_mode::Bool = false
    end

    #=
    function toBSON(params::Param, path::String)
        res = Dict(fn=>getfield(params, fn) for fn ∈ fieldnames(typeof(params)))
        bson(path, res)
    end

    function fromBSON(path::String)::Param
        res = BSON.load(path)
        return Param(
            res[:st_arr],
            res[:en_arr],
            res[:mu_f],
            res[:sigma_f],
            res[:unif_log_lik],
            res[:L],
            res[:max_beta],
            res[:seed],
            res[:n_win],
            res[:mode_arr],
            res[:step_size],
            res[:fixed_inference_flag],
            res[:predef_beta_arr],
            res[:single_end_mode],
            res[:n_max_ats],
            res[:n_min_ats],
            res[:min_ws],
            res[:max_unif_ws],
        )
    end
    =#

    @with_kw mutable struct EMData
        ws::Union{AbstractArray, Nothing} = nothing
        alpha_arr::Union{Vector, Nothing} = nothing
        beta_arr::Union{Vector, Nothing} = nothing
        ats_arr::Union{Vector, Nothing} = nothing
        lb_arr::Union{Vector, Nothing} = nothing
        bic::Union{Number, Nothing} = NaN
        label = nothing
        fragment_size::String = "normal"
    end

    function toDict(self::EMData)::OrderedDict
        res = OrderedDict(
            "ws" => ".",
            "alpha_arr" => ".",
            "beta_arr" => ".",
            "ats_arr" => ".",
        )

        for i = keys(res)
            j =  getfield(self, Symbol(i))
            if isnothing(j)
                j = "NA"
            else
                j = join(map(string, j), ",")
            end

            if j == ""
                continue
            end
            res[i] = j
        end

        res["fragment_size"] = self.fragment_size
        res["bic"] = isnothing(self.bic) ? "NA" : string(self.bic)

        return res
    end

    function EMHeader()::Vector{String}
        temp = toDict(emptyData())
        key = collect(keys(temp))
        return ["utr", key...]
    end

    function emptyData()::EMData
        return EMData()
    end

    struct SplitData
        n_win::Number 
        st_arr::Vector
        en_arr::Vector
        ws_arr::Vector
        mode_arr::Vector
    end

    function dnorm(weights::Union{Number, Vector}; mean::Number, sd::Number, do_log::Bool = false)::Union{Number, Vector}
        if do_log
            return logpdf(Normal(mean, sd), weights)
        end
        return pdf(Normal(mean, sd), weights)
    end

    function cal_z_k(alpha_arr::Vector, beta_arr::Vector, ws::Vector, k::Int, log_zmat; params::Param)::AbstractArray

        K = length(ws) - 1
        
        if 0 < k <= K
            if params.cage_mode
                log_zmat[:, k] = log(ws[k]) .+ lik_l_ab(
                    params.st_arr, alpha_arr[k], beta_arr[k]; 
                    do_log = true
                )
            else
                log_zmat[:, k] = log(ws[k]) .+ lik_lr_ab(
                    params.st_arr, params.en_arr, alpha_arr[k], beta_arr[k]; 
                    params=params, do_log = true
                )
            end
        else
            log_zmat[:, K + 1] .= log(ws[K + 1]) + params.unif_log_lik
        end
        return log_zmat
    end

    # Z is log likelihood
    function norm_z(Z::AbstractArray)::AbstractArray
        tmp_max_vec = maximum(Z, dims=2)
        Z = exp.(Z .- tmp_max_vec)
        return Z ./ sum(Z, dims=2)
    end

    # maximize ws given Z
    function maximize_ws(Z, max_unif_ws::Number)
        ws = sum(Z, dims = 1) / size(Z)[1]

        m = length(ws)
        if ws[m] > max_unif_ws
            ws[1:(m - 1)] = (1 - max_unif_ws) * ws[1:(m - 1)]./ sum(ws[1:(m-1)])
            ws[m] = max_unif_ws
        end
        # 2d row and 1d column error
        return vec(ws)
    end

    function mstep(alpha_arr::Vector, beta_arr::Vector, ws::AbstractArray, Z::AbstractArray, k::Int; params::Param)::EMData
        K = length(ws) - 1 # last component is uniform component
        new_ws = maximize_ws(Z, params.max_unif_ws)

        new_alpha_arr = alpha_arr
        new_alpha_arr[k] = sum(Z[:,k] .* params.st_arr) / sum(Z[:,k])
        
        new_beta_arr = beta_arr;
        temp_beta = sqrt.( sum(Z[:,k] .* (params.st_arr .- new_alpha_arr[k]) .^ 2 ) / sum(Z[:,k]) )

        if isnan(temp_beta)
            temp_beta = params.max_beta
        elseif temp_beta < params.step_size
            temp_beta = params.step_size
        end

        temp_abs = abs.(params.predef_beta_arr .- temp_beta)

        new_beta_arr[k] = params.predef_beta_arr[findall(x -> x == minimum(temp_abs), temp_abs)][1]
    
        if new_beta_arr[k] > params.max_beta
            new_beta_arr[k] = params.max_beta
        end

        return EMData(
            ws = new_ws, 
            alpha_arr = new_alpha_arr, 
            beta_arr = new_beta_arr
        )
    end

    function elbo(log_zmat::AbstractArray, Z::AbstractArray)::Number
        LZ = deepcopy(Z)
        LZ[findall(!isequal(0), Z)] .= log.(Z[findall(!isequal(0), Z)])
        entropy = -1 .* Z .* LZ

        lb = exp_log_lik(log_zmat, Z) .+ sum(entropy)
        if isnan(lb)
            throw(string("lower bounder is na: ", lb))
        end
        return lb
    end

    function exp_log_lik(log_zmat::AbstractArray, Z::AbstractArray)::Number
        ZZ = Z .* log_zmat
        ZZ[findall(isequal(0), Z)] .= 0
        return sum(ZZ)
    end

    function lik_f0(L::Number; do_log::Bool=true)::Number
        if do_log
            return -2 * log(L)
        end
        return 1 / L / L
    end

    function lik_f0_single(mu_f::Number, sigma_f::Number, L::Number; do_log::Bool = true)::Number
        if do_log
            return -log(L) .+ dnorm(mu_f; mean=mu_f, sd=sigma_f, do_log=do_log)
        else
            return 1/L .* dnorm(mu_f; mean=mu_f, sd=sigma_f, do_log=do_log)
        end
    end

    # p(l,r|alpha,beta)
    function lik_lr_ab(l_arr::Vector, r_arr::Vector, alpha::Number, beta::Number; params::Param, do_log::Bool=false)
        if do_log
            return lik_l_ab(l_arr, alpha, beta, do_log=true) .+ lik_r_l(l_arr, r_arr; mu_f=params.mu_f, sigma_f=params.sigma_f, do_log=true)
        else
            return ws .* lik_l_ab(l_arr,alpha, beta) .* lik_r_l(l_arr, r_arr; mu_f=params.mu_f, sigma_f=params.sigma_f)
        end
    end

    # p(r|l)
    function lik_r_l(l_arr::Vector, r_arr::Vector; mu_f::Number, sigma_f::Number, do_log::Bool=false)::Vector
        return dnorm(l_arr .- r_arr; mean = mu_f, sd = sigma_f, do_log = do_log)
    end

    # p(l|alpha, beta)
    function lik_l_ab(l_arr::Vector, alpha::Number, beta::Number; do_log::Bool=false)::Vector
        return dnorm(l_arr; mean = alpha, sd = beta, do_log = do_log)
    end

    function lik_r_s(r_arr::Vector, s::Number, do_log::Bool=false)::Vector
        res = [r_arr .<= s] ./ s
        if do_log
            return log.(res)
        else
            return res
        end
    end

    # generate random k such that each K is a group and no consecutive elements are the same
    function gen_k_arr(K::Int64, n::Int64, seed::Int=2)::Vector{Int64}
        if K == 0
            return zeros(n)
        elseif K == 1
            return [K for _ in 1:n]
        elseif K == 2
            return [i % 2 + 1 for i = 0:(n - 1)]
        else
            nn = ceil(Int, n / K)
            res = zeros(nn * K)
            res[1:K] = StatsBase.sample(MersenneTwister(seed), 1:K, K, replace = false)

            for i = 2:nn
                st = (i - 1) * K
                res[(st + 1):(st + K)] = StatsBase.sample(MersenneTwister(seed), 1:K, K, replace = false)

                if res[st] == res[st + 1]
                    tmpind = StatsBase.sample(MersenneTwister(seed), 2:K, 1, replace = false)[1]
                    tmp = res[st + tmpind]
                    res[st .+ tmpind] = res[st .+ 1]
                    res[st .+ 1] = tmp
                end
            end

            return res[1:n]
        end
    end

    function cal_bic(log_zmat::AbstractArray, Z::AbstractArray)::Number
        N, K = size(Z)
        K = K - 1
        return -2 * exp_log_lik(log_zmat, Z) + (3 * K + 1) * log(N)
    end

    # perform inference given alpha_arr and beta_arr
    function fixed_inference(alpha_arr::Vector, beta_arr::Vector, n_frag::Number, nround::Int; params::Param)::EMData

        lb, lb_arr = -Inf, Array{Number}(undef, nround, 1)
        N = n_frag
        K = length(alpha_arr)

        ws = rand(MersenneTwister(params.seed), K+1) .+ 1
        if K == 0
            ws[1] = ws[1] + 1
        else
            ws[1:K] = ws[1:K] .+ 1
        end
        ws = ws ./ sum(ws)

        k_arr = gen_k_arr(length(alpha_arr), nround, params.seed)
        log_zmat = zeros( N, K + 1)

        for k = 1:K
            log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k, log_zmat; params=params)
        end

        Z = nothing
        i = 1
        for i = 1:nround
            debug(LOGGER, Formatting.format(
                FormatExpr("fixed_inference: iteration={}, lb={}"), i, lb
            ))

            log_zmat = cal_z_k(alpha_arr, beta_arr, ws,k_arr[i], log_zmat, params=params)
            
            Z = norm_z(log_zmat)
            ws = maximize_ws(Z, params.max_unif_ws)

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
        debug(LOGGER, Formatting.format(
            FormatExpr("fixed_inference: Run all {} iterations. lb={}"), 
            i, lb
        ))
        bic = NaN
        if !isnothing(Z)
            bic = cal_bic(log_zmat, Z)
        end
        
        debug(LOGGER, Formatting.format(
            FormatExpr("fixed_inference: bic={}; estimated ws: {}; estimated alpha:  {}; estimated beta: {}"), 
            bic, ws, alpha_arr, beta_arr
        ))

        lb_arr = lb_arr[findall(x -> !isnan(x), lb_arr)]
        return EMData(
            ws = ws, 
            alpha_arr = alpha_arr, 
            beta_arr = beta_arr, 
            lb_arr = lb_arr, 
            bic = bic
        )
    end

    # perform inference for K components
    # theta_arr => theta_win_mat
    function em_algo(ws::Vector, para_mat::EMData, nround::Int; params::Param)::EMData
        lb = -Inf
        lb_arr = [NaN for _ = 1:nround]
        K = length(ws) - 1

        # assume this is a NamedTuple
        alpha_arr = para_mat.alpha_arr
        beta_arr = para_mat.beta_arr

        k_arr = gen_k_arr(length(alpha_arr), nround, params.seed)
        log_zmat = zeros(length(params.st_arr), K+1)

        for k = 1:(K+1)
            log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k, log_zmat; params=params)
        end

        Z = nothing
        i = 1
        
        for i = 1:nround
            log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k_arr[i], log_zmat; params=params)
            Z = norm_z(log_zmat)
            res = mstep(alpha_arr, beta_arr, ws, Z, k_arr[i]; params=params) # , log_zmat
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
        debug(LOGGER, Formatting.format(
            FormatExpr("em_algo: Run all {} interactions. lb={}"), i, lb
        ))

        bic = NaN
        if !isnothing(Z) 
            bic = cal_bic(log_zmat, Z)
        end
        
        debug(LOGGER, Formatting.format(
            FormatExpr("em_algo: bic={}; estimated ws: {}; estimated alpha:  {}; estimated beta: {}"), 
            bic, ws, alpha_arr, beta_arr
        ))

        lb_arr = lb_arr[findall(x -> !isnan(x), lb_arr)]

        sorted_inds = sortperm(alpha_arr)
        alpha_arr = alpha_arr[sorted_inds]
        beta_arr = beta_arr[sorted_inds]
        ws[1:K] = ws[sorted_inds]

        return EMData(
            ws = ws, 
            alpha_arr = alpha_arr, 
            beta_arr = beta_arr, 
            lb_arr = lb_arr, 
            bic = bic
        )
    end

    function rle(value::Vector)
        if length(value) == 0
            return [], []
        end
        k, v = [value[1]], [1]
        for i = 2:length(value)
            
            if value[i] == value[i - 1]
                v[length(v)] += 1
            else
                push!(k, value[i])
                push!(v, 1)
            end
        end
        return v, k
    end

    # calculate coverage 
    function get_border(temp_lens::Vector)
    
        chng = cumsum(temp_lens)
        tmpn = length(chng)

        st_arr = [1, chng[1:(length(chng) - 1)] .+ 1...]
        en_arr = chng

        return st_arr, en_arr
    end

    # handle zero in sign array, half changed to previous/next segment non-zero sign
    function handle_zero(sign_arr::Vector)::Vector
        if sum(findall(x -> x == 0, sign_arr)) == 0
            return sign_arr
        end

        # first sign is always 1
        if sign_arr[1] != 1
            throw(string("sign_arr[1] === 1"))
        end

        temp_lens, temp_vals = rle(sign_arr)
        st_arr, en_arr = get_border(temp_lens)

        zero_inds = findall(x -> x == 0, temp_vals)

        for i = zero_inds
            st = st_arr[i]
            en = en_arr[i]

            if i == length(temp_vals) || en - st <= 1
                sign_arr[st:en] .= temp_vals[i - 1]
            else
                mid = round(Int, (st + en) / 2)
                sign_arr[st:mid] .= temp_vals[i - 1]
                sign_arr[(mid + 1):en] .= temp_vals[i + 1]
            end
        end

        return sign_arr
    end

    function split_data(l_arr::Vector, step_size::Number; npoints::Int = 1024)::SplitData
        # step_size defined in outer function
        # set boudary of kernel density, to keep the last st_arr .<= en_arr
        boundary = (max(1, minimum(l_arr) - 10) + 1, maximum(l_arr) - 1)
        if boundary[2] - boundary[1] < 1024
            npoints = round(Int, float(boundary[2] - boundary[1])) - 1
        end

        ks_res = kde(
            l_arr, bandwidth=5 * step_size, 
            npoints=npoints, kernel = Normal,
            boundary=boundary
        )
        
        # change points
        sign_arr = sign.(diff([-1, ks_res.density...]))
        sign_arr = handle_zero(sign_arr)
        
        lens = rle(sign_arr)
        chng = round.(Int, cumsum( lens[1] ))

        mode_arr = collect(ks_res.x)[ chng[collect(1:2:length(chng))] ]
        n_mode = length(mode_arr)
        boarder_arr = collect(ks_res.x)[ chng[collect(2:2:length(chng))] ][1:(n_mode-1)]

        st_arr = [max(1, minimum(l_arr) - 10), boarder_arr .+ 1...]
        en_arr = [boarder_arr..., maximum(l_arr) + 10]

        ws_arr = zeros(n_mode) .+ 1

        for i = 1:n_mode
            ws_arr[i] = sum(st_arr[i] .<= l_arr .& l_arr .<= en_arr[i])
        end
        ws_arr = ws_arr / sum(ws_arr)

        return SplitData(n_mode, st_arr, en_arr, ws_arr, mode_arr)
    end

    function init_para(n_ats::Int, st_arr::Vector, params::Param)::EMData
        tmp_res = split_data(st_arr, params.step_size)

        if tmp_res.n_win >= n_ats
            tmpinds = StatsBase.sample(MersenneTwister(params.seed), 1:tmp_res.n_win, n_ats, replace = false)
        else
            tmpinds = [
                collect(1:tmp_res.n_win)..., 
                StatsBase.sample(MersenneTwister(params.seed), 1:tmp_res.n_win, Weights(tmp_res.ws_arr), n_ats - tmp_res.n_win, replace=true)...
            ]
        end

        tmpinds = sort(tmpinds)

        res = emptyData()
        res.alpha_arr = zeros(n_ats)
        for i = 1:n_ats
            ti = tmpinds[i]
            tmpst = tmp_res.st_arr[ti]
            tmpen = tmp_res.en_arr[ti]

            if tmpst > tmpen
                throw("tmpst must < tmpen")
            elseif tmpst == tmpen
                res.alpha_arr[i] = tmpst
            else
                res.alpha_arr[i] = rand(
                    MersenneTwister(params.seed), 
                    Uniform(tmpst, tmpen), 1
                 )[1]
            end
        end

        res.beta_arr = rand(MersenneTwister(params.seed), Uniform(params.step_size, params.max_beta), n_ats)

        return res
    end

    function em_optim0(n_ats::Int64; params::Param, n_trial::Int=20)::EMData

        lb_arr = [-Inf for _ = 1:n_trial]
        bic_arr = [-Inf for _ = 1:n_trial]
        res_list = Array{EMData}(undef, n_trial, 1)

        for i = 1:n_trial
            try
                ws = rand(MersenneTwister(params.seed), n_ats + 1) .+ 1
                ws[n_ats + 1] = params.max_unif_ws
                ws = ws / sum(ws)

                # initilize alpha_arr and beta_arr, considered pre_alpha_arr and pre_beta_arr
                para_mat = init_para(n_ats, params.st_arr, params)

                res_list[i] = em_algo(ws, para_mat, n_trial; params=params)
            catch e
                # println(e)
                debug(LOGGER, Formatting.format(
                    FormatExpr("em_optim0: Error found  in {} trial. Next - {}"), i, e
                ))
                res_list[i] = emptyData()
            end

            if isnothing(res_list[i].alpha_arr)
                continue
            end

            lb_arr[i] =  res_list[i].lb_arr[length(res_list[i].lb_arr)]
            bic_arr[i] =  res_list[i].lb_arr[length(res_list[i].bic)]

            debug(LOGGER, Formatting.format(
                FormatExpr("em_optim0: K = {}, {}_trial, n_trial_{}: ws: {}; alpha: {}; beta: {};lb = {}; bic = {}"), 
                n_ats, i, n_trial, round.(res_list[i].ws),
                res_list[i].alpha_arr, res_list[i].beta_arr, 
                res_list[i].lb_arr[length(res_list[i].lb_arr)],
                res_list[i].bic 
            ))
            
        end

        min_ind = argmin(bic_arr)
        return res_list[min_ind]
    end

    function rm_component(res::EMData; n_frag::Number, params::Param, nround::Int=200)

        K = length(res.alpha_arr)
        res.alpha_arr = round.(Int, res.alpha_arr)
        res.beta_arr = round.(Int, res.beta_arr)

        rm_inds = findall(x -> x < params.min_ws, res.ws[1:K])

        if sum(rm_inds) == 0
            return res
        end

        alpha_arr = Vector()
        beta_arr = Vector()

        for i = 1:length(res.alpha_arr)
            if i in rm_inds
                continue
            end
            push!(alpha_arr, res.alpha_arr[i])
            push!(beta_arr, res.beta_arr[i])
        end

        return fixed_inference(alpha_arr, beta_arr, n_frag, nround, params=params)
    end

    function atsmix(params::Param, error_log=nothing)::EMData
        nround = 50
        try
            n_frag = length(params.st_arr)

            if !params.cage_mode && n_frag != length(params.en_arr)
                throw("the length of st_arr != en_arr")
            end
    
            nround=50
    
            lb_arr = [-Inf for _ = 1:params.n_max_ats]
            bic_arr = [Inf for _ = 1:params.n_max_ats]
            res = nothing

            for i = params.n_max_ats:-1:params.n_min_ats
                temp = em_optim0(i, params=params)
    
                if isnothing(res)
                    res = temp
                elseif !isnothing(temp.bic) && !isnothing(res.bic) && !isnan(temp.bic) && !isinf(temp.bic)
                    if temp.bic < res.bic
                        res = temp
                    end
                end
            end
            # println(res)
            debug(LOGGER, string(res))
            if isnothing(res.ws) 
                debug(LOGGER, "fit: Inference failed. No results available.")
                return res
            end
            
            # remove low weight component
            res = rm_component(res, n_frag=n_frag, params=params)
            
            # calculate read assignment
            if isnothing(res.label)
                N = n_frag
                K = length(res.ws) - 1
                log_zmat = zeros(N, K + 1)
                for k = 1:(K + 1)
                    log_zmat = cal_z_k(res.alpha_arr, res.beta_arr, res.ws, k, log_zmat; params=params)
                end
                # Z = norm_z(log_zmat)
                label = argmax.(log_zmat)
                res.label = label
            end
            
            if params.single_end_mode
                res.beta_arr = res.beta_arr .+ params.sigma_f 
            end

            res.fragment_size = "normal"
            len_arr = abs.(params.en_arr .- params.st_arr) .+ 1
            mu_len = mean(len_arr)
            if abs(mu_len - params.mu_f) > 2 * params.sigma_f
                if mu_len - params.mu_f < 0
                    res.fragment_size = "short"
                else
                    res.fragment_size = "long"
                end
            end   
            
            return res
        catch e
            if isnothing(error_log)
                error_log = "error.log"
            end

            open(error_log, "a+") do w
                write(w, string(
                    string(e), "\t",
                    join(map(string, params.st_arr), ","), "\t",
                    join(map(string, params.en_arr), ","), "\t",
                    params.L, "\n"
                ))
                close(w)
            end
            # println(e)
            debug(LOGGER, string(e))
        end
        return emptyData()
    end

    function atsmix_by_R(params::Param, error_log=nothing)::EMData

        script = joinpath(@__DIR__, "atsmix.R")
        len_arr = abs.([y - x + 1 for (x, y) = zip(params.st_arr, params.en_arr)])
        
        n_max_ats = params.n_max_ats
        n_min_ats = params.n_min_ats
        st_arr = params.st_arr
        en_arr = params.en_arr
        L = params.L
        mu_f = params.mu_f
        sigma_f = params.sigma_f
        min_ws = params.min_ws
        max_beta = params.max_beta
        fixed_inference_flag = params.fixed_inference_flag
        single_end_mode = params.single_end_mode
        max_unif_ws = params.max_unif_ws
        seed = params.seed
        cage_mode = params.cage_mode

        try
            atsmix = R"""
            source($script)
            atsmix(n_max_ats=$n_max_ats, 
                    n_min_ats=$n_min_ats,
                    st_arr=$st_arr, 
                    en_arr=$en_arr, 
                    L=$L, 
                    mu_f=$mu_f,
                    sigma_f=$sigma_f,
                    min_ws = $min_ws,
                    max_beta = $max_beta,
                    fixed_inference_flag = $fixed_inference_flag,
                    single_end_mode = $single_end_mode,
                    max_unif_ws = $max_unif_ws,
                    seed = $seed,
                    cage_mode = $cage_mode
            )
            """

            res = rcopy(atsmix)
            res = collect(res)
            ws = res[1][2]
            if isa(ws, Number)
                ws = [ws]
            end

            alpha_arr = res[2][2]
            if isa(alpha_arr, Number)
                alpha_arr = [alpha_arr]
            end

            beta_arr = res[3][2]
            if isa(beta_arr, Number)
                beta_arr = [beta_arr]
            end

            lb_arr = res[4][2]
            if isa(lb_arr, Number)
                lb_arr = [lb_arr]
            end

            label = nothing
            if length(res) >= 6
                label = res[6][2]
            end

            fragment_size = nothing
            if length(res) >= 7
                fragment_size = res[7][2]
            end

            if isnothing(alpha_arr)
                return emptyData()
            end

            return EMData(
                ws = ws, 
                alpha_arr = alpha_arr, 
                beta_arr = beta_arr,
                lb_arr = lb_arr,
                bic = res[5][2],
                label = label,
                fragment_size = fragment_size
            )
        catch e
            
            if isnothing(error_log)
                error_log = "error.log"
            end

            open(error_log, "a+") do w
                write(w, string(
                    string(e), "\t",
                    join(map(string, st_arr), ","), "\t",
                    join(map(string, en_arr), ","), "\t",
                    L, "\n"
                ))
                close(w)
            end
            # println(e)
            debug(LOGGER, string(e))
            return emptyData()
        end
    end

    function fit(
        utr,
        # maximum number of ATS sites
        n_max_ats::Int, 
        # minimum number of ATS sites
        n_min_ats::Int, 

        # information for all DNA fragments, must be n_frag x 1 vector if specified
        # l, start location of each DNA fragment, from 5' to 3'on the 5'-UTR
        st_arr::Vector, 
        # r, end location of each DNA fragment
        en_arr::Vector; 

        # fragment size information
        # fragment length mean
        mu_f::Int = 300, 
        # fragment length standard deviation
        sigma_f::Int = 50, 

        # pa site information
        # minimum weight of ATS site
        min_ws::Float64 = 0.01, 
        # maximum weight of uniform component
        max_unif_ws::Number = 0.1,
        # maximum std for ATS site
        max_beta::Float64 = 50.0,
        seed::Int=42,

        step_size::Number=5,
        # inference with fixed parameters
        fixed_inference_flag::Bool = false,

        #single end mode
        single_end_mode::Bool = false,

        using_R::Bool = true,

        cage_mode::Bool = false,

        exon_coord::Union{Dict{Int, Int}, Nothing} = nothing,
        error_log = nothing,
        debug = false
    )::Union{OrderedDict, EMData}

        if !cage_mode && !single_end_mode && length(st_arr) == length(en_arr)
            len_arr = abs.([y - x + 1 for (x, y) = zip(st_arr, en_arr)])
            mu_f = mean(len_arr)
            sigma_f = std(len_arr)
        end

        L = abs(utr.End - utr.Start)

        unif_log_lik = 0
        if single_end_mode
            st_arr = en_arr .+ mu_f
            unif_log_lik = lik_f0_single(mu_f, sigma_f, L, do_log = true)
        else
            unif_log_lik = lik_f0(L, do_log = true)
        end

        params = Param(
            abs.(st_arr), 
            abs.(en_arr), 
            mu_f, 
            sigma_f, 
            unif_log_lik,
            L, 
            max_beta,
            seed,
            0,
            Vector(),
            step_size,
            fixed_inference_flag,
            collect(step_size:step_size:max_beta),
            single_end_mode, 
            n_max_ats,
            n_min_ats,
            min_ws,
            max_unif_ws,
            cage_mode
        )

        # toBSON(params, "test_R/test.bson")
        runner = using_R ? atsmix_by_R : atsmix

        res = runner(params, error_log)

        if isnothing(res.bic) || isinf(res.bic)
            return ""
        end

        res.ats_arr = []
        if !isnothing(exon_coord) && !isnothing(res.alpha_arr)
            exon_coord = Dict(abs(y) =>x for (x, y) = exon_coord)

            for i = res.alpha_arr
                if haskey(exon_coord, i)
                    push!(res.ats_arr, exon_coord[i])
                else
                    push!(res.ats_arr, utr.Strand == "+" ? utr.Start + abs(i) : utr.End - abs(i))
                end
            end
        elseif !isnothing(res.alpha_arr)
            for i = res.alpha_arr
                push!(res.ats_arr, utr.Strand == "+" ? utr.Start + abs(i) : utr.End - abs(i))
            end
        end

        if debug
            return res
        end

        res = toDict(res)
        res["utr"] = Genomic.get_bed_short(utr)
        return res        
    end
end
