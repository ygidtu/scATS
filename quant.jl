#!/usr/bin/env julia
#=
Created at 2021.01.18

This script is used to quantify the expression counts of each peaks
=#

using ArgParse
using BioAlignments
using FilePathsBase
using GenomicFeatures
using CodecZlib
using Memento
using ProgressMeter
using XAM

#=
using StatsBase
=#
include(joinpath(@__DIR__, "src", "bam.jl"))
include(joinpath(@__DIR__, "src", "genomic.jl"))

logger = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")

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


function load_bed(input_file::String)::Dict{String, Vector{Vector}}
    beds = Dict()

    header = nothing
    open(input_file, "r") do r
        seekend(r)
        fileSize = position(r)

        seekstart(r)
        p = Progress(fileSize, 1)   # minimum update interval: 1 second
        while !eof(r)
            update!(p, position(r))
            line = split(strip(readline(r)), "\t")

            if isnothing(header)
                header = line
                continue
            end
            
            temp = Dict(x => y for (x, y) in zip(header, line))

            if occursin(",",  temp["gene_name"])
                update!(p, position(r))
                continue
            end

            site = split(temp["utr"], r"[:-]")

            strand = site[length(site)]
            if strand != "+"
                strand = "-"
            end

            chrom, start_pos, end_pos = site[1], parse(Int, site[2]), parse(Int, site[3])
            alpha = split(get(temp, "inferred_site", get(temp, "infered_sites", "")), ",")
  
            if length(alpha) > 0
                try
                    if !haskey(beds, temp["gene_name"])
                        beds[temp["gene_name"]] = Vector()
                    end

                    sites = Vector()

                    for x = alpha
                        if x != ""
                            s = round(Int, parse(Float64, x))
                            
                            if strand == "+"
                                push!(sites, Genomic.BED(chrom, s, s + 1, temp["gene_name"], string(length(beds) + 1), strand))
                            else
                                push!(sites, Genomic.BED(chrom, s - 1, s, temp["gene_name"], string(length(beds) + 1), strand))
                            end
                        end
                    end

                    push!(beds[temp["gene_name"]], sites)
                catch e
                    error(logger, e)
                end
            end

            update!(p, position(r))
        end
        close(r)
    end

    return beds
end


function generate_region_by_vec_of_regions(regions::Vector)
    regions = sort(regions)

    return Genomic.BED(
        regions[1].Chrom, 
        regions[1].Start - 1, 
        regions[length(regions)].End + 1, 
        "", "", regions[1].Strand
    )
end


function main(input_file::String, bam::String, output::String)
    beds = load_bed(input_file)

    if length(beds) < 1
        info(logger, "There is no region needs to counts")
        exit(0)
    end

    bam_list = Bam.prepare_bam_list(bam)

    if length(bam_list) < 1
        info(logger, "There is no bam file")
        exit(0)
    end

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

    l = ReentrantLock()

    # Counting by multi-threads
    p = Progress(length(beds), 1, "Counting...")

    wc = open(string(output, ".count.gz"), "w+")
    wp = open(string(output, ".psi.gz"), "w+")

    sc = GzipCompressorStream(wc)
    sp = GzipCompressorStream(wp)

    Threads.@threads for name = collect(keys(beds))
        utrs = beds[name]

        sums = Dict{String, Int}()

        counts = Vector()
        
        # counting
        for bs = utrs
            res = Bam.count_reads(
                bam_list,
                generate_region_by_vec_of_regions(bs)
            )

            push!(counts, res)

            lock(l)
            for b = bs
                bStr = Genomic.get_bed_key(b)
                site = b.Start
                if b.Strand != "+"
                    site = b.End
                end

                for (key, value) in get(res, site, Dict())
                    if value > 0
                        sums[key] = get(sums, key, 0) + value
                        write(sc, string(bStr, "\t", key, "\t", value, "\n"))
                    end
                end
            end
            flush(wc)
            unlock(l)
        end

        # PSI
        for i in 1:length(utrs)
            bs =  utrs[i]
            res = counts[i]

            for b = bs
                bStr = Genomic.get_bed_key(b)
                site = b.Start
                if b.Strand != "+"
                    site = b.End
                end

                lock(l)
                for (key, value) in get(res, site, Dict())
                    if value > 0
                        write(sp, string(bStr, "\t", key, "\t", value / sums[key], "\n"))
                    end
                end
                flush(wp)
                unlock(l)
            end
        end
        
        next!(p)
    end

    # force to write string to gzip file
    # to avoid unexcept end of gzip file
    write(sc, "\n")
    write(sp, "\n")

    close(sc)
    close(sp)
    close(wc)
    close(wp)
end


ccall(:jl_exit_on_sigint, Nothing, (Cint,), 0)
# main(args["input"], args["bam"], args["output"])
try
    main(args["input"], args["bam"], args["output"])
catch ex
    println(ex)
    if isa(ex, InterruptException)
        info(logger, "caught keyboard interrupt")
    end
end
