using Random
using StatsBase
using Memento
using Compat: @__MODULE__  # requires a minimum of Compat 0.26. Not required on Julia 0.7

# Create our module level logger (this will get precompiled)
const LOGGER = getlogger(@__MODULE__)

# Register the module level logger at runtime so that folks can access the logger via `get_logger(MyModule)`
# NOTE: If this line is not included then the precompiled `MyModule.LOGGER` won't be registered at runtime.
function __init__()
    Memento.register(LOGGER)
end


Number = Union{Int64, Float64}

function replace_with_closest(ref_arr, query_arr)
    n = length(query_arr)
    res = zeros(n)

    for i = 1:n
        tmp_ind = argmin(abs(ref_arr .- query_arr[i]))
        res[i] = ref_arr[tmp_ind]
    end
    return res
end


function gen_k_arr(K::Int64, n::Int64)::Vector{Int64}

    if K <= 1
        return [K - 1 for _ in 1:n]
    else
        res = Vector()
        count_index = 0
        pool = [x for x in 0:(K - 1)]
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


function rle_wrap(x::Vector)
    vals, lens = rle(x)

    index = Vector()
    for i = vals
        push!(index, findfirst(isequal(i), x))
    end
    return index, lens, vals
end