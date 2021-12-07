module ATSModel
    using Distributions
    using Formatting
    using KernelDensity
    using Memento
    using OrderedCollections
    using Parameters
    using Peaks
    using Random
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

    function setLevel(level::String)
        setlevel!(LOGGER, level)
    end

    Number = Union{Int64, Float64, Int}

    export init, setLevel, getHeader, toDict, run

    ## Parameters and related functions
    @with_kw mutable struct Param
        st_arr::Vector
        L::Int = 500   # Length of UTR
        n_max_ats::Int = 5
        n_min_ats::Int = 1
        min_ws::Number = 0.01
        max_unif_ws::Number = 0.1
        max_beta::Number = 50
        fixed_inference_flag::Bool = false
        seed::Int = 42

        # initialize and other info
        step_size::Int = 5
        nround::Int = 50
        unif_log_lik::Union{Nothing, Number} = nothing
        predef_beta_arr::Union{Nothing, Vector{Number}} = nothing

        # Result
        utr = nothing # Genomic.Region
        gene_name::String = ""
        transcript_name::String = ""
        num_of_reads::Int = 0 # number of reads
        inferred_sites::Vector{Int} = []
        alpha_arr::Vector{Number} = []
        beta_arr::Vector{Number} = []
        ws::Vector{Number} = []
        bic::Float64 = -Inf
        lb_arr::Vector{Number} = []
    end

    function init(
        st_arr::Vector,
        utr;
        n_max_ats::Int=5, 
        n_min_ats::Int=1, 
        min_ws::Float64=0.01, 
        max_unif_ws::Float64 = 0.1, 
        fixed_inference_flag::Bool = false,
        step_size::Int = 5,
        nround::Int = 50,
        max_beta::Number = 50,
        seed::Int = 42
    )::Param
        return Param(
            st_arr = st_arr, 
            L = utr.End - utr.Start,
            utr = utr,
            num_of_reads = length(st_arr),
            gene_name = utr.Name,
            transcript_name = utr.Score, 
            n_max_ats = n_max_ats, 
            n_min_ats =  n_min_ats,
            min_ws = min_ws, 
            max_unif_ws = max_unif_ws, 
            max_beta = max_beta,
            fixed_inference_flag = fixed_inference_flag,
            step_size = step_size, 
            nround = nround,
            seed = seed
        )
    end

    # Results and related functions
    function toDict(self::Union{Param, Nothing})::OrderedDict
        res = OrderedDict(
            "utr" => nothing,
            "gene_name" => ".",
            "transcript_name" => ".",
            "num_of_reads" => 0,
            "inferred_sites" => ".",
            "alpha_arr" => ".",
            "beta_arr" => ".",
            "ws" => ".",
            "L" => 0,
            # "bic" => 0
        )

        if isnothing(self)
            return res
        end

        site = []
        if self.utr.Strand == "+"
            site = self.utr.Start .+ self.alpha_arr
        else
            site = self.utr.End .- self.alpha_arr
        end

        return OrderedDict(
            "utr" => string(self.utr.Chrom, ":", self.utr.Start, "-", self.utr.End, ":", self.utr.Strand),
            "gene_name" => self.gene_name,
            "transcript_name" => self.transcript_name,
            "num_of_reads" => self.num_of_reads,
            "inferred_sites" => join(site, ","),
            "alpha_arr" => join(self.alpha_arr, ","),
            "beta_arr" => join(self.beta_arr, ","),
            "ws" => join(self.ws, ","),
            "L" => self.L,
            "bic" => self.bic
        )
        return res
    end

    function getHeader()::Vector{String}
        temp = toDict(nothing)
        return collect(keys(temp))
    end

    # Model related functions

    ## Evaluate the probability density
    function lik_l_ab(l_arr::Vector, alpha::Number, beta::Number, do_log::Bool = false)::Vector{Number}
        norm = Normal(alpha, beta)
        if do_log
            return logpdf.(norm, l_arr)
        end
        return pdf.(norm, l_arr)
    end

    function cal_z_k(self::Param, k::Int, log_zmat)
        if 0 < k <= length(self.alpha_arr)
            log_zmat[:, k] = log.(self.ws[k]) .+ lik_l_ab(
                self.st_arr, 
                self.alpha_arr[k],
                self.beta_arr[k],
                true
            )
        else
            log_zmat[:, k] = log.(self.ws[k]) .+ self.unif_log_lik
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
    function maximize_ws(self::Param, Z)::Vector{Number}
        ws = sum(Z, dims = 1) / size(Z)[1]

        m = length(ws)
        if ws[m] > self.max_unif_ws
            ws[1:(m - 1)] = (1 - self.max_unif_ws) * ws[1:(m - 1)]./ sum(ws[1:(m-1)])
            ws[m] = self.max_unif_ws
        end
        # 2d row and 1d column error
        return vec(ws)
    end

    function mstep(self::Param, Z::AbstractArray, k::Int)
        tmp_sumk = sum(Z[:,k])

        if tmp_sumk < 1e-8
            Z[:,k] = Z[:,k] .+ 1e-8
            Z = norm_z(Z)
            tmp_sumk = sum(Z[:,k])
        end

        self.ws = maximize_ws(self, Z)
        return self
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
            return -1 * log(L)
        end
        return 1 / L
    end

    # generate random k such that each K is a group and no consecutive elements are the same
    function gen_k_arr(K::Int64, n::Int64)::Vector{Int64}

        if K <= 1
            return ones(n)  
            # Python version is zeros here, but this is index of array, therefore change to ones now
        end

        ii = 0
        last_ind = -1
        arr = randperm(K)
        res = Vector()
        for i = 1:n
            if ii % K == 0
                if arr[1] == last_ind
                    tmpi = rand(1:K)
                    arr[1], arr[tmpi] = arr[tmpi], arr[1]
                end

                ii = 0
                last_ind = arr[length(arr)]
            end
            push!(res, arr[ii + 1])
            ii += 1
        end

        return res
    end

    function cal_bic(log_zmat::AbstractArray, Z::AbstractArray)::Number
        N, K = size(Z)
        K = K - 1
        return -2 * exp_log_lik(log_zmat, Z) + (3 * K + 1) * log(N)
    end

    # perform inference given alpha_arr and beta_arr
    ## Param or Results
    function fixed_inference(self::Param)::Param
        self.ws = init_ws(self, length(self.alpha_arr))
        return em_algo(self)
    end

    # perform inference for K components
    # theta_arr => theta_win_mat
    function em_algo(self::Param)::Param
        lb = -Inf
        lb_arr = []
        N = length(self.st_arr)
        K = length(self.alpha_arr) - 1

        k_arr = gen_k_arr(length(self.alpha_arr), self.nround)
        log_zmat = zeros(N, K+1)

        for k = 1:(K+1)
            log_zmat = cal_z_k(self, k, log_zmat)
        end

        Z = nothing
        i = 1
        
        for i = 1:self.nround
            # E-Step
            log_zmat = cal_z_k(self, k_arr[i], log_zmat)
            Z = norm_z(log_zmat)

            self = mstep(self, Z, k_arr[i]) # , log_zmat
            
            lb_new = elbo(log_zmat, Z)
            push!(lb_arr, lb_new)
            lb_arr[i] = lb_new
            
            if isinf(lb_new)
                lb = Inf
                break
            end

            if abs(lb_new - lb) < abs(1e-6 * lb)
                break
            else
                lb = lb_new
            end
        end
        
        self.lb_arr = lb_arr[findall(x -> !isnan(x), lb_arr)]

        sorted_inds = sortperm(self.alpha_arr)
        self.alpha_arr = self.alpha_arr[sorted_inds]
        self.beta_arr = self.beta_arr[sorted_inds]
        self.ws = self.ws[sorted_inds]
        if !isnothing(Z) 
            self.bic = cal_bic(log_zmat, Z)
        end

        return self
    end

    function sample_alpha(self::Param, n_ats::Int)::Vector{Number}
        # 随机生成alpha值，先根据高斯分布随机生成高斯分布数据，然后找峰
        boundary = (max(1, minimum(self.st_arr) - 10) + 1, maximum(self.st_arr) - 1)
        npoints = (self.L + 100) + 100

        kernel = kde(
            self.st_arr,
            npoints=npoints, 
            kernel = Normal,
            boundary=boundary
        )

        x_arr, y_arr = kernel.x, kernel.density

        peak_inds, _ = findmaxima(y_arr)
        
        if length(peak_inds) < 1
            return []
        end

        peaks = x_arr[peak_inds]
        peaks_ws = y_arr[peak_inds] ./ sum(y_arr[peak_inds])

        # 峰太多了就随机挑选，太少了就补
        if n_ats <= length(peaks)
            return StatsBase.sample(MersenneTwister(self.seed), peaks, ProbabilityWeights(peaks_ws), n_ats, replace = false)
        else
            mu = StatsBase.sample(MersenneTwister(self.seed), peaks, ProbabilityWeights(peaks_ws), n_ats - length(peaks), replace = true)

            push!(peaks, mu...)
            shift = round.(Int, rand(MersenneTwister(self.seed), Normal(0, 5 * self.step_size), n_ats)) # MersenneTwister(42), 

            return peaks .+ shift
        end
    end

    function init_ws(self::Param, n_ats::Int)
        ws = rand(Uniform(), n_ats)

        ws = ws ./ sum(ws)

        if ws[length(ws)] > self.max_unif_ws
            ws[1:(length(ws) - 1)] = ws[1:(length(ws) - 1)] .* (1-self.max_unif_ws)
            ws[length(ws)] = self.max_unif_ws
        end
        return ws
    end

    # init alpha_arr, beta_arr, ws
    function init_para(self::Param, n_ats::Int)::Param
        self.alpha_arr = sample_alpha(self, n_ats)
        self.beta_arr = StatsBase.sample(MersenneTwister(self.seed), self.predef_beta_arr, n_ats, replace = true)
        self.ws = init_ws(self, n_ats)
        return self
    end

    # remove sites with weights < minimum weight
    function rm_component(self::Param)::Param
        keep_inds = findall(x -> x >= self.min_ws, self.ws)

        if length(keep_inds) == length(self.ws)
            return self
        end

        self.alpha_arr = self.alpha_arr[keep_inds]
        self.beta_arr = self.beta_arr[keep_inds]
        self.ws = []
        return fixed_inference(self)
    end

    function em_optim0(self::Param, n_ats::Int)::Param
        n_trial = 5

        res_list = []
        for _ = 1:n_trial
            para = init_para(deepcopy(self), n_ats)
            para = em_algo(para)
            push!(res_list, deepcopy(para))
        end

        min_ind = argmin([x.bic for x = res_list])
        return res_list[min_ind]
    end

    function run(self::Param)
        if length(self.st_arr) < 1
            return nothing
        end

        if self.max_beta < self.step_size
            error(LOGGER, Formatting.format(
                FormatExpr("max_beta = {}, step_size = {}, max_beta has to be greater than step_size!"), 
                self.max_beta, self.step_size
            ))
            return nothing
        end

        self.predef_beta_arr = collect(self.step_size:self.step_size:self.max_beta)

        self.unif_log_lik = lik_f0(self.L, do_log = true)
        res_list = []
        for n_ats = self.n_max_ats:-1:self.n_min_ats
            # use deepcopy to save data, or the res will be override
            res = em_optim0(self, n_ats)
            push!(res_list, deepcopy(res))
        end

        # choose the smallest BIC value
        min_ind = argmin([x.bic for x = res_list])
        res = rm_component(res_list[min_ind])
        # convert alpha_arr from float to int
        res.alpha_arr = round.(Int, res.alpha_arr)
        return res
    end

end
