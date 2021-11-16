#!/usr/bin/env julia
#=
Created at 2021.01.18

This script is used to quantify the expression counts of each peaks
=#

using ArgParse
using BioAlignments
using FilePathsBase
using GenomicFeatures
using Libz
using Memento
using ProgressMeter
using XAM

#=
using StatsBase
=#
include(joinpath(@__DIR__, "src", "bam.jl"))
include(joinpath(@__DIR__, "src", "genomic.jl"))


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
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

    return parse_args(s)
end

args = parse_commandline()


function load_bed(input_file::String)
    beds = []  # Dict{String, Vector}()

    header = nothing
    open(input_file, "r") do r
        seekend(r)
        fileSize = position(r)

        seekstart(r)
        p = Progress(fileSize, 1)   # minimum update interval: 1 second
        while !eof(r)

            line = split(strip(readline(r)), "\t")

            if isnothing(header)
                header = line
                continue
            end
            
            temp = Dict(x => y for (x, y) in zip(header, line))

            site = split(temp["utr"], r"[:-]")

            strand = site[length(site)]
            if strand != "+"
                strand = "-"
            end

            chrom, start_pos, end_pos = site[1], parse(Int, site[2]), parse(Int, site[3])
            alpha = split(get(temp, "inferred_site", get(temp, "infered_sites", "")), ",")

            if length(alpha) > 0
                try
                    for x = alpha
                        if x != ""
                            s = round(Int, parse(Float64, x))
                            
                            if strand == "+"
                                push!(beds, Genomic.BED(chrom, s, s + 1, temp["gene_name"], string(length(beds) + 1), strand))
                            else
                                push!(beds, Genomic.BED(chrom, s - 1, s, temp["gene_name"], string(length(beds) + 1), strand))
                            end
                        end
                    end
                catch e
                end
            end

            update!(p, position(r))
        end
        close(r)
    end

    return beds
end


function main(input_file::String, bam::String, output::String)
    beds = load_bed(input_file)

    bam_list = Bam.prepare_bam_list(bam)

    key = basename(dirname(dirname(bam)))
    key = split(string(key), "-")[1]

    # bam = joinpath("/mnt/raid64/Covid19_Gravida/apamix/bam", string(key, ".bam"))

    output = absolute(Path(output))
    out_dir = parent(output)
    try
        if !exists(out_dir)
            mkdir(Path(out_dir), recursive=true)
        end
    catch e
        error(logger, Formatting.format(FormatExpr("Error while create {}: {}"), output, e))
        exit(1)
    end
    output = string(output)

    open(output, "w+") do w
        stream = nothing

        if endswith(output, ".gz")
            stream = ZlibDeflateOutputStream(w)
        end

        p = Progress(length(beds), 1, "Computing...")
        Threads.@threads for i = eachindex(beds)
            b = beds[i]
            res = Bam.reads(bam_list, b)

            for (key, value) in res
                if value > 0
                    write(isnothing(stream) ? w : stream, string(Genomic.get_bed_key(b), "\t", key, "\t", value, "\n"))
                end
            end
            next!(p)

            b = nothing
            res = nothing
            
            # GC.safepoint()

            # if i % 100 ==  0
            #     GC.gc()
            # end
        end

        if !isnothing(stream)
            close(stream)
        end

        close(w)
    end
end


main(args["input"], args["bam"], args["output"])
