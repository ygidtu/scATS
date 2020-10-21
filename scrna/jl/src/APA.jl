
module APA
    using Distributions
    using Formatting
    using KernelDensity
    using StatsBase


    include("EM.jl")
    include("utils.jl")

    mutable struct APAData
        n_max_apa::Int64
        n_min_apa::Int64
        r1_utr_st_arr::Vector{Number}
        r1_len_arr::Vector{Number}
        r2_len_arr::Vector{Number}
        polya_len_arr::Vector{Number}
        pa_site_arr::Vector{Number}
        utr_len::Vector{Number}
        LA_dis_arr::Vector{Number}
        pmf_LA_dis_arr::Vector{Number}
        pmf_s_dis_arr::Vector{Number}
        min_LA::Int64
        max_LA::Int64
        mu_f::Int64
        sigma_f::Int64
        min_ws::Float64
        min_pa_gap::Int64
        max_beta::Int64
        pre_alpha_arr::Vector{Number}
        pre_beta_arr::Vector{Number}
        theta_step::Int64
        cb::String
        umi::String
        fixed_inference_flag::Bool
        pre_init_theta_win_mat
        region_name::String
        verbose::Bool
        all_theta::Vector{Number}
        n_all_theta::Int64
        n_frag::Int64
        nround::Int64
        lb_arr::Vector{Number}
        bic_arr::Vector{Number}
        ws_flag_arr::Vector{Number}
        res_lst::Vector
        all_theta_lik_mat
        all_theta_lik_mat_list::Vector
        pre_init_alpha_arr::Vector{Number}
        init_beta_arr::Vector{Number}
        unif_log_lik::Vector{Number}
    end

    mutable struct APAWins
        st_arr::Vector{Number}
        en_arr::Vector{Number}
        ws_arr::Vector{Number}
        mode_arr::Vector{Number}
    end

    mutable struct InitThetaWin
        res_alpha_arr::Vector{Number}
        res_beta_arr::Vector{Number}
    end

    Base.show(io::Io, self::InitThetaWin) print(
        io,
        "init alpha: ", string(self.res_alpha_arr),
        "\ninti beta: ", string(self.res_beta_arr)
    )

    function check(self::APAData)::APAData
        if all(isnan.(self.pa_site_arr)) && length(self.pa_site_err) != length(self.r1_len_arr)
            self.pa_site_err = [NaN for _ in 1:length(self.r1_len_arr)]
        end

        if all(isnan.(self.polya_len_arr)) && length(self.polya_site_err) != length(self.r1_len_arr)
            self.polya_site_err = [NaN for _ in 1:length(self.r1_len_arr)]
        end

        n_pre_pa_sites = 0
        if length(self.pre_alpha_arr) > 0
            self.min_theta = min(minmum(self.r1_len_arr), minmum(self.pre_alpha_arr))
            self.pre_alpha_arr = sort(self.pre_alpha_arr)
            if maximum(self.pre_alpha_arr) > self.L
                throw("The max of pre_alpoha_arr was out of UTR length")
            end
            n_pre_pa_sites = length(self.pre_alpha_arr)
            self.pre_init_alpha_arr = replace_with_closest(self.all_theta, self.pre_alpha_arr)
        else
            self.pre_init_alpha_arr = Vector()
        end

        if length(self.pre_beta_arr) > 0
            if n_pre_pa_sites != length(self.pre_beta_arr)
                throw(string("The num of pre_pa_sites (", n_pre_pa_sites, ") wasn't same with pre_beta_arr (", length(self.pre_beta_arr), ")."))
            end
            
            # if pre_beta_arr is all NaN
            if sum(isnan.(self.pre_beta_arr)) == length(self.pre_beta_arr)
                throw("All pre_beta_arr was NaN.")
            end

            if maximum(self.pre_beta_arr) > self.max_beta
                throw(string("The max of pre_beta_arr (", maximum(self.pre_beta_arr), ") was more than ", self.max_beta))
            end
        end

        if n_pre_pa_sites > 0
            if n_pre_pa_sites > self.n_max_apa
                warn(LOGGER, string("n_max_apa: {", self.n_max_apa, "} change to n_pre_pa_sites: {", n_pre_pa_sites, "}"))
                self.n_max_apa = n_pre_pa_sites
                self.lb_arr = [-Inf for _ = 1:self.n_max_apa]
                self.bic_arr = [-Inf for _ = 1:self.n_max_apa]
                self.ws_flag_arr = [false for _ = 1:self.n_max_apa]                
            end

            if n_pre_pa_sites < self.n_min_apa
                wawrn(LOGGER, string("n_min_apa: {", self.n_max_apa, "} change to n_pre_pa_sites: {", n_pre_pa_sites, "}"))
                self.n_min_apa = n_pre_pa_sites
            end
        end

        if self.max_beta < self.theta_step
            throw(string("max_beta: {", self.max_beta, "} wasn't more than {", self.theta_step, "}"))
        end

        if self.n_frag != length(self.r1_utr_st_arr) || self.n_frag != length(self.r2_len_arr)
            throw(string("n_frag (", self.n_frag, ") wasn't same to r1_utr_st_arr (", self.r1_utr_st_arr, ")"))
        end

        if self.L < self.min_theta
            throw(string("L (", self.L, ") wasn't more than min_theta (", self.min_theta, ")"))
        end
        return self
    end

    function new(
        n_max_apa::Int64,
        n_min_apa::Int64,
        r1_utr_st_arr::Vector{Number},
        r1_len_arr::Vector{Number},
        r2_len_arr::Vector{Number},
        polya_len_arr::Vector{Number},
        pa_site_arr::Vector{Number},
        utr_len::Vector{Number},
        LA_dis_arr::Vector{Number}=[10, 30, 50, 70, 90, 110, 130],
        pmf_LA_dis_arr::Vector{Number}=[309912, 4107929, 802856, 518229, 188316, 263208, 101],
        min_LA::Int64=20,
        max_LA::Int64=150,
        mu_f::Int64=300,
        sigma_f::Int64=50,
        min_ws::Float64=.01,
        min_pa_gap::Int64=100,
        max_beta::Int64=70,
        pre_alpha_arr::Vector{Number}=Vector(),
        pre_beta_arr::Vector{Number}=Vector(),
        theta_step::Int64=9,
        cb::String="",
        umi::String="",
        fixed_inference_flag::Bool=false,
        pre_init_theta_win_mat=nothing,
        region_name::String="",
        verbose::Bool=false
    )::APAData
        all_theta = [x for x in min_theta:theta_step:L]

        if length(LA_dis_arr) < 1
            s_dis_arr = [x for x in min_LA:10:max_LA]
            pmf_LA_dis_arr = [1 / length(s_dis_arr) for _ in 1:length(s_dis_arr)]
        end

        pmf_s_dis_arr = pmf_LA_dis_arr ./ sum(pmf_LA_dis_arr)

        self = APAData(
            n_max_apa,
            n_min_apa,
            r1_utr_st_arr,
            r1_len_arr,
            r2_len_arr,
            polya_len_arr,
            pa_site_arr,
            utr_len,
            LA_dis_arr,
            pmf_LA_dis_arr,
            pmf_s_dis_arr,
            min_LA,
            max_LA,
            mu_f,
            sigma_f,
            min_ws,
            min_pa_gap,
            max_beta,
            pre_alpha_arr,
            pre_beta_arr,
            theta_step,
            cb,
            umi,
            fixed_inference_flag,
            pre_init_theta_win_mat,
            region_name,
            verbose,
            all_theta, 
            length(all_theta),
            length(r1_len_arr), 
            50,
            [-Inf for _ = 1:n_max_apa],
            [Inf for _ = 1:n_max_apa],
            [false for _ = 1:n_max_apa],
            Vector(), 
            zeros(length(r1_len_arr), 
            length(all_theta)), 
            Vector(), 
            Vector(),
            Vector(),
            EM.lik_f0(utr_len, max_LA, true)
        )

        return init(check(self))
    end

    function newRes()::EM.EMData
        return EM.EMData(NaN, Vector(), Vector(), Vector(), NaN)
    end

    function kernal_smooth(self::APAData, beta::Vector{Number})::Vector{Number}
        n_step = beta / self.theta_step
        n_bd_step = ceil(n_step * 3)
        weights = pdf(Normal(0, n_step), [i for i = -n_bd_step:(n_bd_step) - 1])
        weights = weights ./ sum(weights)

        n_all_theta = length(self.all_theta)

        if  n_all_theta != size(self.all_theta_lik_mat)[2]
            throw(string("n_all_theta: {", n_all_theta, "} wasn't equal to column of all_theta_lik_mat"))
        end

        all_win_log_lik_mat = zeros(self.n_frag, n_all_theta)

        if sum(self.all_theta_lik_mat) > 0
            for i = 1:n_bd_step
                tmp_inds = [x for x in max(i - n_bd_step, 1):(i + n_bd_step)]
      
                tmp_weight = weights[tmp_inds[tmp_inds .> 0]]
                tmp_weight = tmp_weight ./ sum(tmp_weight)
 
                try
                    all_win_log_lik_mat[:, i] = log(sum(tmp_weight * self.all_theta_lik_mat[:, tmp_inds], dims=2))
                catch e
                    error(logger, string(e))
                    exit(1)
                end
            end

            for i = n_bd_step:(n_all_theta - n_bd_step)
                tmp_inds = [x for x = (i - n_bd_step):(i + n_bd_step)]
                tmp_weight = weight
                all_win_log_lik_mat[:, i] = log(sum(tmp_weight * self.all_theta_lik_mat[:, tmp_inds], dims=2))
            end

            for i = (n_all_theta - n_bd_step):n_all_theta
                tmp_inds = [x for x = (i - n_bd_step):(i + n_bd_step)]
                tmp_weight = weights[tmp_inds < n_all_theta]
                tmp_weight = tmp_weight ./ sum(tmp_weight)

                tmp_inds = tmp_inds[tmp_inds < n_all_theta]

                all_win_log_lik_mat[:, i] = log(sum(tmp_weight * self.all_theta_lik_mat[:, tmp_inds], dims=2))
            end
        end
        return all_win_log_lik_mat
    end

    function split_data(self::APAData)::APAWins
        n_frag = length(self.r1_utr_st_arr)
        coverage_cnt = zeros(L)

        for i = 1:n_frag
            if self.r1_utr_st_arr[i] >= 1
                tmp_inds = self.r1_utr_st_arr[i] - 1 + np.arange(self.r1_len_arr[i])
                coverage_cnt[tmp_inds] += 1
            end
        end

        ks_res_x = [x for x in 1:self.L]
        ks_res = kde((coverage_cnt, ks_res_x))
        ks_res = ks_res.y

        sign_arr = sign.(diff(append!([-1], ks_res)))

        if sign_arr[1] != 1
            throw("First sign value not 1, split_data func")
        end

        if any(sign_arr .== 0)
            st_arr, lengths, tmp_vals = rle_wrap(sign_arr) 
            en_arr = cumsum(lengths)
            tmp_n = length(tmp_vals)

            for (index, value) in enumerate(tmp_vals)
                if value == 0
                    st = st_arr[index]
                    en = en_arr[index]
                    if  index == tmp_n - 1 || en - st <= 1
                        sign_arr[st:en] .= tmp_vals[i - 1]
                    else
                        mid = round((st + en + 1) / 2)
                        sign_arr[st:mid] = tmp_vals[i - 1]
                        sign_arr[mid:en] = tmp_vals[i + 1]
                    end
                end
            end
        end
        st_arr, lengths, tmp_vals = rle_wrap(sign_arr)
        chng = cumsum(lengths)
        mode_arr = ks_res_x[chng[x for x in 1:2:length(chng)] - 1]
        n_mode = length(mode_arr)

        boarder_arr = ks_res_x[chng[x for x in 1:2:length(chng)] - 1][[x for x in 1:n_mode]]

        st_arr = append!([1], boarder_arr .+ 1)
        en_arr = append(boarder_arr, [L])
        ws_arr = ones(n_mode)

        for i = 1:n_mode
            ws_arr[i] = sum(coverage_cnt[st_arr[i]:en_arr[i]])
        end
        ws_arr = ws_arr  ./ sum(ws_arr)

        return APAWins(
            st_arr, en_arr, ws_arr, mode_arr
        )

    
        function init(self::APAData)
            oth_inds = isnan.(self.pa_site_arr)
            pa_inds = 1 .- oth_inds

            for i = 1:self.n_all_theta
                self.all_theta_lik_mat[oth_inds, i] = EM.lik_lsr_t(
                    self.r1_utr_st_arr[oth_inds],
                    self.r1_len_arr[oth_inds],
                    self.r2_len_arr[oth_inds],
                    self.polya_len_arr[oth_inds],
                    self.all_theta[i],
                    self.mu_f,
                    self.sigma_f,
                    self.pmf_s_dis_arr,
                    self.s_dis_arr
                )
            end

            if sum(pa_inds) > 0
                for i in findall(isequal(1), pa_inds)
                    tmp_ind = argmin(abs(self.all_theta .- self.pa_site_arr[i]))
                    self.all_theta_lik_mat[i, tmp_ind] = EM.lik_lsr_t0(
                        self.r1_utr_st_arr[i],
                        self.r1_len_arr[i],
                        self.r2_len_arr[i],
                        self.polya_len_arr[i],
                        self.all_theta[tmp_ind],
                        self.pmf_s_dis_arr,
                        self.s_dis_arr,
                        false
                    )
                end
            end

            self.init_beta_arr = [x for x in self.theta_step:self.theta_step:ceil(self.max_beta / self.theta_step) * self.theta_step + 1]

            for i = 1:length(self.init_beta_arr)
                push!(self.all_win_log_lik_mat_list, kernal_smooth(self, self.init_beta_arr[i]))
            end
            return self
        end
    end

    function sample_theta(k_arr::Vector{Number}, all_theta::Vector{Number}, data_wins::APAWins, mode::String="full")::Vector{Number}
        res_arr = [.0 for _ in 1:length(k_arr)]
        for (i, iw) = enumerate(k_arr)
            if mode == "full"
                tmp_theta_arr = all_theta[
                    findall(
                        x -> data_wins.st_arr[iw] <= x < data_wins.en_arr[iw],
                        all_theta
                    )
                ]
            elseif mode == "left"
                tmp_theta_arr = all_theta[
                    findall(
                        x -> data_wins.st_arr[iw] <= x < data_wins.mode_arr[iw],
                        all_theta
                    )
                ]
            elseif mode == "right"
                tmp_theta_arr = all_theta[
                    findall(
                        x -> data_wins.mode_arr[iw] <= x < data_wins.en_arr[iw],
                        all_theta
                    )
                ]
            else
                throw(string("Unknown mode: ", mode))
            end

            if length(tmp_theta_arr) > 0
                res_arr[i] = sample(tmp_theta_arr, 1, replace = false)
            else
                res_arr[i] = -1
            end
        end
        return sort(res_arr)
    end

    function sample_theta1(k_arr::Vector{Number}, n_data_win::Number, all_theta::Vector{Number}, data_wins::APAWins)::Vector{Number}
        tmp_ind = k_arr .== n_data_win
        theta_arr1 = Vector()
        if sum(tmp_ind) > 0
            theta_arr1 = sample_theta([n_data_win], all_theta, data_wins, "full")
        end

        tmp_ind = k_arr .!= n_data_win
        theta_arr2 = Vector()
        if sum(tmp_ind) > 0
            theta_arr2 = sample_theta([n_data_win], all_theta, data_wins, "left")
        end
        append!(theta_arr1, theta_arr2)
        return theta_arr1
    end

    function get_data_win_lab(alpha_arr::Vector{Number}, data_wins::APAWins)::Vector{Number}
        n_data_win = length(data_wins.st_arr)
        left = 1
        right = data_wins.mode_arr[1]
        lab_arr = copy(alpha_arr)
        for i = 1:n_data_win
            lab_arr[findall(x -> left <= x < right, alpha_arr)] = i
            left = data_wins.mode_arr[i]

            right = nothing
            if i < n_data_win - 1
                right = data_wins.mode_arr[i+1]
            end
        end
        lab_arr[alpha_arr .>= left] = n_data_win
        return lab_arr
    end

    function init_theta_win(self::APAData, n_apa::Int64, n_max_trial::Int64=200)::InitThetaWin
        if length(self.pre_alpha_arr) > 0
            pre_init_alpha_arr = replace_with_closest(
                self.all_theta, self.pre_alpha_arr
            )
            pre_init_beta_arr = sample(
                self.init_beta_arr[[x for x = 1:round(Int, length(self.init_beta_arr) / 2 + 1)]], 
                length(self.pre_alpha_arr), replace = true
            )

            if length(self.pre_beta_arr) > 0
                tmp_inds = 1 .- isnan.(self.pre_alpha_arr)
                pre_init_beta_arr[tmp_inds] = replace_with_closest(
                    self.init_beta_arr, self.pre_beta_arr[tmp_inds]
                )
            end

            if length(self.pre_alpha_arr) == n_apa
                return InitThetaWin(pre_init_alpha_arr, pre_init_beta_arr)
            end
        end

        k_big_inds = findfirst(x -> x > 0.1, self.data_wins.ws_arr)
        n_data_win = length(self.data_wins.ws_arr)
        k1 = length(k_big_inds)
        
        tmp_alpha_arr = Vector{Number}()
        for i = 1:n_max_trial
            if i == n_max_trial
                throw(string("Failed to generate valid theta after {", n_max_trial, "} trial."))
            end

            if n_data_win >= n_apa
                k_arr = sort(sample(n_data_win, n_apa, self.data_wins.ws_arr, replace = false))

                tmp_alpha_arr = sample_theta(k_arr, self.all_theta, self.data_wins, "right")
            elseif n_data_win >= n_apa - k1
                theta_arr1 = sample_theta([x for x = 1:n_data_win], self.all_theta, self.data_wins, "right")

                if k1 == 1
                    k_arr = k_big_inds .+ 1
                else
                    ws_tmp = self.data_wins.ws_arr[k_big_inds]
                    k_arr = sample(k_big_inds .+ 1, n_apa - n_data_win, ws_tmp ./ ws_tmp.sum(), replace = false)
                end

                theta_arr2 = sample_theta1(k_arr, n_data_win, self.all_theta, self.data_wins)

                tmp_alpha_arr = sort(append!(theta_arr1, theta_arr2))
            else
                theta_arr1 = sample_theta([x for x = 1:n_data_win], self.all_theta, self.data_wins, "right")
                theta_arr2 = sample_theta1(k_big_inds .+ 1, n_data_win, self.all_theta, self.data_wins)

                k_arr = sample(n_data_win, n_apa - n_data_win - k1, self.data_wins.ws_arr)

                theta_arr3 = sample_theta(k_arr, self.all_theta, self.data_wins, "full")

                append!(theta_arr1, theta_arr2)
                append!(theta_arr1, theta_arr3)
                tmp_alpha_arr = sort(theta_arr1)
            end
            flag1, flag2 = all(diff(tmp_alpha_arr .>= self.min_pa_gap)), all(tmp_alpha_arr .> 0)

            if flag1 && flag2
                break
            end
        end

        if length(self.pre_alpha_arr) > 0
            sam_lab = get_data_win_lab(tmp_alpha_arr, self.data_wins)
            pre_lab = get_data_win_lab(self.pre_alpha_arr, self.data_wins)

            n_to_rm = 0
            for i = 1:length(pre_lab)
                smp_inds = findall(x -> 0 < x && x == pre_lab[i], sam_lab)
                
                if sum(smp_inds) == 0
                    n_to_rm += 1
                elseif sum(smp_inds) == 1
                    sam_lab[smp_inds] = 0
                else
                    tmp_ind = argmin(abs.(tmp_alpha_arr[smp_inds] - self.pre_alpha_arr[i]))
                    sam_lab[smp_inds[tmp_ind]] = 0
                end
            end
            if n_to_rm > 0
                smp_inds = sam_lab ./ 0
                if sum(smp_inds) > 0
                    throw("spm_inds out of range, continue next trial")
                end
                tmp_inds = sample(smp_inds, n_to_rm)
                sam_lab[tmp_inds] = 0
            end

            valid_inds = sam_lab .> 0

            res_alpha_arr = copy(tmp_alpha_arr[valid_inds])
            append!(res_alpha_arr, self.pre_init_alpha_arr)
       
            return InitThetaWin(
                sort(res_alpha_arr),
                sample(
                    self.init_beta_arr[[x for x = 1:round(Int, length(self.init_beta_arr) / 2 + 1)]], 
                    n_apa
                )
            )
        
        return InitThetaWin(
            tmp_alpha_arr,
            sample(
                self.init_beta_arr[[x for x = 1:round(Int, length(self.init_beta_arr) / 2 + 1)]], 
                n_apa
            )
        )
    end

    function em_optim(self::APAData, n_apa::Int64, n_trial::Int64=5, verbose::Bool=false)::EM.EMData
        lb_arr = [-Inf for i in 1:n_trial]
        bic_arr = [Inf for i in 1:n_trial]

        res_list = [nothing for _ = 1:n_trial]

        for i = 1:n_trial
            try
                ws = [rand(Uniform()) + 1 for _ in 1:n_apa]
                ws[n_apa] += -1
                ws = ws ./ sum(ws)
                theta_win_mat = init_theta_win(self, n_apa)

                if verbose
                    debug(logger, string(theta_win_mat))
                end

                res_list[i] = EM.em_algo(
                    ws,
                    theta_win_mat,
                    self.n_frag,
                    self.pre_alpha_arr,
                    self.pre_beta_arr,
                    self.all_theta,
                    self.init_beta_arr,
                    self.all_win_log_lik_mat_list,
                    self.pre_init_alpha_arr,
                    self.min_pa_gap,
                    self.unif_log_lik,
                    self.min_theta,
                    self.L,
                    verbose
                )
            catch e
                debug(logger, string("Error found  in ", i, "trial. Next"))
                res_list[i] = newRes()
            end

            if any(diff(res_list[i].alpha_arr .< self.min_pa_gap))
                continue
            end

            lb_arr[i] =  res_list[i].lb_arr[length(res_list[i].lb_arr)]
            bic_arr[i] =  res_list[i].lb_arr[length(res_list[i].bic)]

            if verbose
                debug(logger, Formatting.format(
                    FormatExpr("K = {}, {}_trial, n_trial_{}: ws: {}\nalpha: {}\nbeta: {}\nlb = {}\nbic = {}"), 
                    n_apa, i, n_trial, round.(res_list[i].ws),
                    res_list[i].alpha_arr, res_list[i].beta_arr, 
                    res_list[i].lb_arr[length(res_list[i].lb_arr)],
                    res_list[i].bic 
                ))
            end
        end

        min_ind = argmin(bic_arr)
        return res_list[min_ind]
    end

    function rm_component(
        res::EM.EMData, min_ws::Int64, 
        pre_alpha_arr::Vector{Number},
        pre_init_alpha_arr::Vector{Number},
        n_frag::Int,
        all_theta::Vector{Number},
        all_win_log_lik_mat_list::Vector{Number},
        init_beta_arr::Vector{Number},
        unif_log_lik,
        nround::Int=200,
        verbose::Bool=false)::EM.EMData

        K = length(res.alpha_arr)
        if K < 1
            return res
        end

        kept_inds = Vector()
        if length(pre_alpha_arr) > 0
            pre_flag = [x in pre_init_alpha_arr for x = res.alpha_arr]
            kept_inds = findall(isequal(0), [x < min_ws && 1 - y for (x, y) in zip(res.ws[1:K], pre_flag)])
        else
            kept_inds = findall(isqeual(1), res.ws[1:K] .>= min_ws)
        end

        if length(kept_inds) == K || length(kept_inds) == 0
            # log_mat = zeros(n_frag,  K + 1)

            # for k in 1:(K+1)
            #     loz_zmat = EM.cal_z_k(
            #         all_win_log_lik_mat_list,
            #         res.alpha_arr,
            #         res.beta_arr,
            #         res.ws,
            #         k,
            #         log_zmat,
            #         all_theta,
            #         init_beta_arr,
            #         unif_log_lik
            #     )
            # end
            # z = EM.norm_z(log_zmat)
            # label = argmax(log_zmat, dims=2)
            # res.label = label
            return res
        end

        alpha_arr = res.alpha_arr[kept_inds]
        beta_arr = res.beta_arr[kept_inds]

        return EM.fixed_inference(
            alpha_arr,
            beta_arr,
            n_frag,
            nround,
            all_theta,
            all_win_log_lik_mat_list,
            init_beta_arr,
            unif_log_lik,
            verbose
        )
    end


    function run(self::APAData)::
        data_wins = split_data(self)

        if self.fixed_inference_flag && !isnothing(self.pre_init_theta_win_mat)

            alpha_arr = replace_with_closest(
                self.all_theta, 
                self.pre_init_theta_win_mat[:, 1]
            )

            beta_arr = replace_with_closest(
                self.all_theta, 
                self.pre_init_theta_win_mat[:, 2]
            )

            res = EM.fixed_inference(
                alpha_arr, self.n_frag,
                self.nround, self.all_theta,
                self.all_win_log_lik_mat_list,
                self.init_beta_arr, self.unif_log_lik,
                self.verbose
            )

            res.alpha_arr = self.pre_init_theta_win_mat[:, 1]
            res.beta_arr = self.pre_init_theta_win_mat[:, 2]

            return res
        end

        res_list = Vector{EM.EMData}()
        lb_arr = [-Inf for i in 1:self.n_max_apa]
        bic_arr = [Inf for i in 1:self.n_max_apa]

        for i in self.n_max_apa:-1:self.n_min_apa
            res_tmp = em_optim(self, i, self.verbose)
            push!(res_list, res_tmp)

            if length(res_tmp.ws) < 1
                continue
            end

            lb_arr[i] = res_tmp[lb_arr][length(res_tmp[lb_arr])]
            bic_arr[i] = res_tmp.bic
        end

        res_list = reverse(res_list)
        min_ind = argmin(bic_arr)
        res = res_list[mid_ind]

        if length(res) < 1
            warn(logger, "Inference failed. No results available.")
            return newRes()
        end

        return rm_component(
            res, 
            self.min_ws, 
            self.pre_alpha_arr, 
            self.pre_init_alpha_arr, 
            self.n_frag, 
            self.all_theta, 
            self.all_theta_lik_mat_list, 
            self.init_beta_arr, 
            self.unif_log_lik, 
            verbose
        )
    end
end

using CSV
using DataFrames

df = DataFrame(CSV.File(
    "/mnt/raid61/Personal_data/zhangyiming/code/afe/tests/1_6274198_6276648.txt",
    delim='\t', header = 0
))

rename(df, ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"])


offset = 150
df.V1 = df.V1 .+ offset
df.V7 = df.V7 .+ offset
L = maximum(df.V1) + maximum(df.V2) + 50
# logger.debug('init running')
test = APA.new(
    5,1,df.V1,df.V2,
    df.V3,df.V6,df.V7,
    L,true
)

res = APA.run(test)

println(res)
# # logger.debug('init infer')
# a = test.inference()
