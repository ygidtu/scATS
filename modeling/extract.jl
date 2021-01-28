#!/usr/bin/env julia

using ArgParse
using BSON
using Distributed


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Path to utr bed file"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Prefix of output file"
            arg_type = String
            required = true
        "--process", "-p"
            help = "How many processes to use"
            arg_type = Int64
            default = 1
        "bams"
            help = "Path to bam"
            required = true
            action = :store_arg
            nargs = '+'
    end

    return parse_args(s)
end

args = parse_commandline()

if args["process"] > 1
    addprocs(args["process"])
end

@everywhere begin
    using FilePathsBase
    using Memento
    using ProgressMeter
    

    include(joinpath(@__DIR__, "src", "extract.jl"))
    include(joinpath(@__DIR__, "src", "genomic.jl"))
end


function main(input_file::String, bams::Vector, output::String)
    logger = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")

    beds = []  # Dict{String, Vector}()
    r = open(input_file, "r")
    seekend(r)
    fileSize = position(r)
    seekstart(r)
    p = Progress(fileSize, 1)   # minimum update interval: 1 second
    while !eof(r)
        temp_bed = Genomic.new_bed(readline(r)) #, chromosomes=chromosomes)
        push!(beds, temp_bed)
        update!(p, position(r))
    end
    close(r)

    # output = absolute(Path(output))
    # try
    #     if !exists(output)
    #         mkdir(Path(output), recursive=true)
    #     end
    # catch e
    #     error(logger, Formatting.format(FormatExpr("Error while create {}: {}"), output, e))
    #     exit(1)
    # end

    info(logger, string("The number of utrs: ", length(beds)))
    res = @showprogress "Computing... " pmap(beds) do b
    # for b = beds
        data = Extract.get_record_from_bam(
            bams, b, 
            cage_mode=false, 
            single_end_mode=false
        )
        return Extract.toDict(data)
    end

    data = Dict()
    for row = res
        for (k, v) = row
            data[k] = v
        end
    end

    bson(string(output), data)
end

main(args["input"], args["bams"], args["output"])

