using CSV
using DataFrames
using BenchmarkTools

dfs = [CSV.read("testfile/0", DataFrame, header = false) for i in readdir("testfile")]

function entropy(df)
    df = df ./ ncol(df)

    ent = 0
    for row in eachrow(df)
        row = values(row)
        row = row .* log.(row)

        ent -= sum(filter(!isnan, row))
    end

    return ent
end


@btime [entropy(x) for x = dfs]

