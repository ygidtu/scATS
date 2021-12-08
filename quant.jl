#!/usr/bin/env julia
#=
Created at 2021.01.18

This script is used to quantify the expression counts of each peaks
=#

using CodecZlib
using FilePathsBase
using ProgressMeter


const jobs = RemoteChannel(()->Channel{Tuple}(32));
const results = RemoteChannel(()->Channel{Tuple}(32));


@everywhere begin
    include(joinpath(@__DIR__, "src", "bam.jl"))
    include(joinpath(@__DIR__, "src", "genomic.jl"))

    function format_bStr(b, site::Int)::String
        return string(
            b.Chrom, ":", 
            ifelse(b.Strand == "+", site, site - 1),
            "-",
            ifelse(b.Strand == "+", site + 1, site),
            ":", b.Strand
        )
    end

    function do_work(jobs, results)
        while true
            utrs, bam_list = take!(jobs)

            sums = Dict{String, Int}()
            counts = Vector()
            cStrs, pStrs = [], []

            # counting
            for bs = utrs
                res = Bam.count_reads(bam_list, bs)

                push!(counts, res)

                for site = bs.Sites

                    bStr = format_bStr(bs, site)

                    for (key, value) in get(res, site, Dict())
                        if value > 0
                            sums[key] = get(sums, key, 0) + value
                            push!(cStrs, string(bStr, "\t", key, "\t", value, "\n"))
                        end
                    end
                end
            end

            # PSI
            for i in 1:length(utrs)
                bs =  utrs[i]
                res = counts[i]
 
                for site = bs.Sites
                    bStr = format_bStr(bs, site)

                    for (key, value) in get(res, site, Dict())
                        if value > 0
                            push!(pStrs, string(bStr, "\t", key, "\t", value / sums[key], "\n"))
                        end
                    end
                end
            end

            put!(results, (cStrs, pStrs))
        end
    end
end


## count related functions
function load_bed(input_file::String)::Dict{String, Vector}
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
            sites = split(get(temp, "inferred_sites", get(temp, "infered_sites", "")), ",")

            if length(sites) > 0
                if !haskey(beds, temp["gene_name"])
                    beds[temp["gene_name"]] = Vector()
                end
                
                sites = sort(round.(Int, parse.(Float64, sites)))

                push!(
                    beds[temp["gene_name"]], 
                    Genomic.new(
                        chrom, 
                        sites[1] - 1, 
                        sites[length(sites)] + 1, 
                        strand,
                        sites = sites
                    )
                )
            end

            update!(p, position(r))
        end
        close(r)
    end

    return beds
end


function make_jobs(beds::Dict, bams::Vector)
    for i = values(beds)
        put!(jobs, (i, bams))
    end
end


function quantification(args::Dict; logger)
    input_file, bam, output = args["input"], args["bam"], args["output"]
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

    # Counting by multi-threads
    p = Progress(length(beds), 1, "Counting...")

    wc = open(string(output, ".count.gz"), "w+")
    wp = open(string(output, ".psi.gz"), "w+")

    sc = GzipCompressorStream(wc)
    sp = GzipCompressorStream(wp)

    # multi-processing
    errormonitor(@async make_jobs(beds, bam_list))

    for p = workers()
        remote_do(do_work, p, jobs, results)
    end

    n = length(beds)
    @elapsed while n > 0
        counts, psis = take!(results)
        for c = counts
            write(sc, c)
            flush(sc)
        end

        for p = psis
            write(sp, p)
            flush(sp)
        end
        next!(p)
        n = n - 1
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