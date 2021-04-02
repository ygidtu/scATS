using ArgParse
using CSV
using DataFrames
using Memento
using ProgressMeter
using RCall


function load_group(path::String)
    res = Dict()
    open(path) do r
        while !eof(r)
            line = readline(r)

            line = split(strip(line), "\t")

            res[line[1]] = line[2]
        end
    end

    return res
end


function main(data::String, group::String, output::String)
    logger = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")

    info(logger, "load  quant data")
    data = rcopy(R"""readRDS($data)""") # CSV.read(data, DataFrame)

    info(logger, "load group")
    group = load_group(group)
    # println(group)

    # println(data)

    # first get sum
    rowSums = Dict(x => Dict() for x = keys(group))
    colnames = names(data)
    colnames = colnames[1:length(colnames)]
    colnames = filter(x -> endswith(x, "SJ_out"),colnames)

    p = Progress(size(data)[1], 1, "Computing sums...")
    Threads.@threads for idx in rownumber.(eachrow(data))
        rowname = string(idx) # row.Column1
        next!(p)
        if !haskey(group, rowname)
            continue
        end

        row = data[idx, :]
        
        key = group[rowname]
        if !haskey(rowSums, key)
            rowSums[key] = Dict()
        end
        for colname = colnames
            val = row[Symbol(colname)]
            try
                # println(val)
                # println(get(rowSums[key], colname, 0))
                rowSums[key][colname] = val + get(rowSums[key], colname, 0)
            catch e
                rowSums[key] = Dict()
                rowSums[key][colname] = val + get(rowSums[key], colname, 0)
            end
        end
    end


    # second calcualate percentage

    open(output, "w+") do w
        write(w, string(",",  join(colnames, ","), "\n"))

        p = Progress(size(data)[1], 1, "Computing sums...")
        Threads.@threads for idx in rownumber.(eachrow(data)) # 
            rowname = string(idx) # row.Column1
            next!(p)

            if !haskey(group, rowname)
                continue
            end

            row = data[idx, :]

            key = group[rowname]
            new_row = [string(key, "|", rowname)]

            for colname = colnames
                val = row[Symbol(colname)]
                sumVal = get(rowSums[key], colname, 0)

                push!(new_row, string(sumVal > 0 ? val / sumVal : 0))
            end

            write(w, string(join(new_row, ","), "\n"))
        end
    end

end



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Path to quant matrix"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Path to output file"
            arg_type = String
            required = true
        "--group", "-g"
            help="path to group information from postprocess.jl"
            arg_type = String
            required = true
    end

    return parse_args(s)
end

args = parse_commandline()


#=
main(
    "/mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/CG.quant",
    "/mnt/raid64/ATS/Personal/zhangyiming/stats/overall_ATS_cage_R_gene",
    "/mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/CG.psi"
)

import os
from glob import glob
from subprocess import check_call

g = "/mnt/raid64/ATS/Personal/zhangyiming/stats/overall_ATS_cage_R_gene"
for f in glob("/mnt/raid64/ATS/Personal/zhangyiming/quant/R_sel/*.quant"):
    o = f.replace(".quant", ".psi")

    check_call(f"julia -t 20 /mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/psi.jl -i {f} -g {g} -o {o}", shell=True)
=#

main(args["input"], args["group"], args["output"])
