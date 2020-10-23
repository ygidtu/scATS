#!/usr/bin/env julia
#=
Created by Zhang Yiming at 2020.10.15

I try to convert apamix from python to julia
=#
using ArgParse
using Distributed


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Path to STAR output directory"
            arg_type = String
            required = true
        "--gtf", "-g"
            help = "Path to output file"
            arg_type = String
            required = true
        "--junctions", "-j"
            help = "Path to formatted junction directory"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Path to output file"
            arg_type = String
            required = true
        "--length", "-l"
            help = "Length of sequenced reads"
            arg_type = Int64
            default = 150
        "--process", "-p"
            help = "How many processes to use"
            arg_type = Int64
            default = 1
    end

    return parse_args(s)
end

args = parse_commandline()
# if args["process"] > 1
#     addprocs(args["process"])
# end

@everywhere begin
    using FilePathsBase
    using Formatting
    using ProgressMeter
    include(joinpath(@__DIR__, "start_psi.jl"))
end


function main(input::String, junctions::String, gtf::String, output::String, len::Int64, process::Int64=1)

    output = absolute(Path(output))
    if !exists(output)
        mkdir(output)    
    end

    fs = Vector()
    for i = readdir(input)
        if endswith(i, ".Aligned.sortedByCoord.out.bam")
            push!(fs, StartPSI.new(
                string(joinpath(input, i)),
                string(gtf), 
                string(joinpath(junctions, replace(i, ".Aligned.sortedByCoord.out.bam" => ".SJ.out.tab"))),
                string(joinpath(output, replace(i, ".Aligned.sortedByCoord.out.bam" => ""))),
                len
            ))
        end
    end

    # if process > 1
    @showprogress pmap(fs) do d
        try
            psi = StartPSI.start(d)
        catch e
            println(e)
        end
    end
    # else
    #     @showprogress for d = fs
    #         psi = StartPSI.run(d)
    #     end
    # end
end

main(
    args["input"], args["junctions"],
    args["gtf"], args["output"], 
    args["length"]
)
