#!/usr/bin/env julia
#=
Created at 2021.12.07

This script is used to quantify the expression counts of each peaks
=#

using ArgParse
using FilePathsBase
using Memento
using ProgressMeter

#=
using StatsBase
=#

# subcommand
s = ArgParseSettings()
@add_arg_table s begin
    "cmd"
        help = "ats or count"
        arg_type = String
        required = true
end;

# ats mode
ats_settings = ArgParseSettings()
@add_arg_table ats_settings begin
    "--bam", "-b"
        help = "Path to bam list"
        arg_type = String
        required = true
    "--gtf", "-g"
        help = "Path to reference annotation file, GTF format"
        arg_type = String
        required = true
    "--output", "-o"
        help = "The path to output file"
        arg_type = String
        required = true
    "--utr-length", "-u"
        help = "The length of UTR region"
        arg_type = Int
        default = 500
    "--threads", "-t"
        help = "How many threads to use"
        arg_type = Int64
        default = 1
    "--n-max-ats"
        help = "the maximum of ats inside a utr"
        arg_type = Int64
        default = 5
    "--n-min-ats"
        help = "the minimum of ats inside a utr"
        arg_type = Int64
        default = 1
    "--min-ws"
        help = "min ws"
        arg_type = Float64
        default = 0.1
    "--max-unif-ws"
        help = "maximum uniform ws"
        arg_type = Float64
        default = 0.1
    "--max-beta"
        help = "maximum beta"
        arg_type = Float64
        default = 50.0
    "--step-size"
        help = "step size"
        arg_type = Int
        default = 5
    "--nround"
        help = "number of round to test"
        arg_type = Int
        default = 50
    "--min-reads"
        help = "minimum reads to construct ATS"
        arg_type = Int
        default = 5
    "--seed"
        help = "seed for consistance results"
        arg_type = Int
        default = 42
    "--fixed-inference"
        help="inference with fixed parameters"
        action = :store_true
end;

# count mode
count_settings = ArgParseSettings()
@add_arg_table count_settings begin
    "--input", "-i"
        help = "Path to merged peaks bed"
        arg_type = String
        required = true
    "--bam", "-b"
        help = "Path to bam list"
        arg_type = String
        required = true
    "--output", "-o"
        help = "Prefix of output file"
        arg_type = String
        required = true
end


if ARGS[1] == "ats"
    include(joinpath(@__DIR__, "ats.jl"))
else
    include(joinpath(@__DIR__, "ats.jl"))
end

logger = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")


function main(subcommand::String)
    args = ARGS[2:length(ARGS)]

    if subcommand == "ats"
        infer(parse_args(args, ats_settings), logger = logger)
    elseif subcommand == "count"
        args = parse_args(args, count_settings)
        quantification(args["input"], args["bam"], args["output"], logger = logger)
    else
        info(logger, "please set running subcommand: ats or count")
        println(Base.stderr, usage_string(s))
        exit(1)
    end
end


ccall(:jl_exit_on_sigint, Nothing, (Cint,), 0)

# try
#     main(ARGS[1])
# catch ex
#     if isa(ex, InterruptException)
#         info(logger, "caught keyboard interrupt")
#     else
#         showerror(stdout, e, catch_backtrace())
#     end
# end

main(ARGS[1])