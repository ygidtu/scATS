#!/usr/bin/env julia
#=
Created by Zhang Yiming at 2020.10.15

I try to convert apamix from python to julia
=#
using ArgParse
using Distributed
using FilePathsBase


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Path to STAR output directory"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Path to output file"
            arg_type = String
            required = true
        "--process", "-p"
            help = "How many processes to use"
            arg_type = Int64
            default = 1
    end

    return parse_args(s)
end

args = parse_commandline()
addprocs(args["process"])

@everywhere begin
    using ProgressMeter

    function format_bed(input_file::String, output_file::String)
        count = 1
        open(output_file, "w+") do w
            open(input_file, "r") do io
                while !eof(io)
                    line = split(strip(readline(io)), r"\s")
                    
                    st = string(parse(Int64, line[2]) - 21)
                    en = string(parse(Int64, line[3]) + 20)
    
                    strand = "-"
                    if line[4] == "1"
                        strand = "+"
                    end
    
                    val = ([line[1], st, en, string("JUNCBJ", count), line[7], strand, st, en, "255,0,0", "2", "20,20", "0,300"])
                    
                    write(w, join(val, "\t"))
                    write(w, "\n")
                    count += 1
                end
                close(io)
            end
            close(w)
        end
    end
end

function main(input::String, output::String)
    output = absolute(Path(output))
    if !exists(output)
        mkdir(output)    
    end

    fs = Dict{String, String}()

    for i = readdir(input)
        if endswith(i, "SJ.out.tab")
            fs[joinpath(input, i)] = joinpath(output, i)
        end
    end

    @showprogress pmap(collect(keys(fs))) do i
        format_bed(i, fs[i])
    end
    rmprocs()
end

main(args["input"], args["output"])
