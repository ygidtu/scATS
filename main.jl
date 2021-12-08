#!/usr/bin/env julia
#=
Created at 2021.12.07

This script is used to quantify the expression counts of each peaks
=#

using ArgParse
using Distributed
using FilePathsBase
using Memento


# subcommand
settings = Dict(
    "ats" => ArgParseSettings(),
    "count" => ArgParseSettings()
)
s = ArgParseSettings()
@add_arg_table s begin
    "cmd"
        help = "ats or count"
        arg_type = String
        required = true
end;

# ats mode
@add_arg_table settings["ats"] begin
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
    "--process", "-p"
        help = "Number of processes to used."
        arg_type = Int
        default = 1
end;

# count mode
@add_arg_table settings["count"] begin
    "--input", "-i"
        help = "Path to merged peaks bed."
        arg_type = String
        required = true
    "--bam", "-b"
        help = "Path to bam list."
        arg_type = String
        required = true
    "--output", "-o"
        help = "Prefix of output file."
        arg_type = String
        required = true
    "--process", "-p"
        help = "Number of processes to used."
        arg_type = Int
        default = 1
end

ccall(:jl_exit_on_sigint, Nothing, (Cint,), 0)
logger = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")


if length(ARGS) < 1 || !haskey(settings, ARGS[1])
    info(logger, "please set running subcommand: ats or count")
    println(Base.stderr, usage_string(s))
    exit(1)
end

subcommand, args = ARGS[1], ARGS[2:length(ARGS)]
args = parse_args(args, settings[subcommand])
addprocs(get(args, "process", 1); exeflags=`--project=$(Base.active_project())`)

if subcommand == "ats"
    include(joinpath(@__DIR__, "ats.jl"))
else
    include(joinpath(@__DIR__, "quant.jl"))
end


try
    if subcommand == "ats"
        infer(args, logger = logger)
    elseif subcommand == "count"
        quantification(args, logger = logger)
    end
catch e
    if isa(e, InterruptException)
        info(logger, "caught keyboard interrupt")
    else
        showerror(stdout, e, catch_backtrace())
    end
end

