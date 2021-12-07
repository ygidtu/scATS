#!/usr/bin/env julia
#=
Created at 2021.01.18

This script is used to quantify the expression counts of each peaks
=#

using ArgParse
using BioAlignments
using CodecZlib
using FilePathsBase
using GenomicFeatures
using Memento
using ProgressMeter
using XAM

#=
using StatsBase
=#
include(joinpath(@__DIR__, "src", "bam.jl"))
include(joinpath(@__DIR__, "src", "genomic.jl"))

## count related functions
function load_bed(input_file::String)::Dict{String, Vector{Vector}}
    beds = Dict()

    header = nothing
    open(input_file, "r") do r
        seekend(r)
        fileSize = position(r)

        seekstart(r)
        p = Progress(fileSize, 1, "Loading...")   # minimum update interval: 1 second
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

            chrom, _, _ = site[1], parse(Int, site[2]), parse(Int, site[3])
            alpha = split(get(temp, "inferred_site", get(temp, "infered_sites", "")), ",")
  
            if length(alpha) > 0
                try
                    if !haskey(beds, temp["gene_name"])
                        beds[temp["gene_name"]] = Vector()
                    end

                    alpha = sort(trunc.(Int, parse.(Float64, x)))

                    push!(
                        beds[temp["gene_name"]], 
                        Genomic.Sites(
                            chrom, 
                            alpha[1] - 1, 
                            alpha[length(sites) + 1], 
                            strand,
                            sites
                        )
                    )
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


function quantification(input_file::String, bam::String, output::String; logger)
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
            res = Bam.count_reads(bam_list, bs)

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