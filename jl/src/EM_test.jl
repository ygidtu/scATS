using Test

include("EM.jl")

function test_lik_f0()
    println(EM.lik_f0(40, 50, false))
    println(EM.lik_f0(40, 50, true))
end


function test_lik_r_s()
    r_arr = [0.4, 0.3, .2, .5, .6, .7, .8, .9, .65, .55]
    println(r_arr)
    println(EM.lik_r_s(r_arr, .5))
end


function test_maximize_ws()
    Z = [[0.1, 0.2, 0.3] [0.6, .5, .4] [.8, .7, .9]]
    println(EM.maximize_ws(Z))
    println(EM.norm_z(Z))
end

function test_exp_log_lik()
    Z = [[0.1, 0.2, 0.3] [0.6, .5, .4] [.8, .7, .9]]
    Z = transpose(Z)
    log_zmat = [log.(x) for x in Z]
    println(EM.exp_log_lik(log_zmat, Z))

    println(EM.elbo(log_zmat, Z))
end


# test_lik_f0()
# test_lik_r_s()
test_maximize_ws()
test_exp_log_lik()
