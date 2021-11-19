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

            if occursin(",",  temp["gene_name"])
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

    temp = Path(string(output, ".temp"))
    done = joinpath(temp, "COUNTING_DONE")

    # counts reads
    if !exists(done)
        if exists(temp)
            rm(temp, recursive=true)
        end
    
        mkdir(Path(temp), recursive=true)
        
        p = Progress(length(beds), 1, "Counting...")
        Threads.@threads for b = beds
            res = Bam.count_reads(bam_list, b)

            open(joinpath(temp, b.Name), "a+") do w
                stream = ZlibDeflateOutputStream(w)
                
                for (key, value) in res
                    if value > 0
                        write(stream, string(Genomic.get_bed_key(b), "\t", key, "\t", value, "\n"))
                    end
                end
                
                close(stream)
                close(w)
            end
        
            next!(p)
        end

        # Add Progress log
        open(done, "w+") do w
            write(w, "DONE")
            close(w)
        end
    end

    # merging results and calculate PSI
    fs = readdir(temp)
    p = Progress(length(fs), 1, "PSI...")

    wc = open(string(output, ".count.gz"), "w+")
    wp = open(string(output, ".psi.gz"), "w+")

    sc = ZlibDeflateOutputStream(wc)
    sp = ZlibDeflateOutputStream(wp)

    l = ReentrantLock()
    Threads.@threads for f = fs

        if string(f) == filename(done)
            continue
        end
        
        counts = Vector{String}()
        data = Dict{String, Dict{String, Int}}()
        sums = Dict{String, Int}()

        for line in eachline(open(joinpath(temp, f)) |> ZlibInflateInputStream)
            line = strip(line)
            push!(counts, line)
            line = split(line, "\t")
            if length(line) != 3
                continue
            end
            
            # log the counts of ATS
            if !haskey(data, line[1])
                data[line[1]] = Dict()
            end
            data[line[1]][line[2]] = parse(Int, line[3])

            # log the sum of ATS
            sums[line[2]] = get(sums, line[2], 0) + parse(Int, line[3])
        end

        lock(l)
        for i in counts
            write(sc, string(i, "\n"))
        end

        # save psi
        for (ats, barcodes) in data
            for (bc, val) in barcodes
                write(sp, string(ats, "\t", bc, "\t", val / sums[bc], "\n"))
            end
        end
        unlock(l)
        next!(p)
    end

    close(sc)
    close(sp)
    close(wc)
    close(wp)

    if exists(temp)
        FilePathsBase.rm(temp, recursive=true, force = true)
    end
end

ccall(:jl_exit_on_sigint, Nothing, (Cint,), 0)
# main(args["input"], args["bam"], args["output"])
try
    main(args["input"], args["bam"], args["output"])
catch ex
    println(ex)
    if isa(ex, InterruptException)
        info(LOGGER, "caught keyboard interrupt")
    end
end
