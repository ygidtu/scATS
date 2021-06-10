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
include(joinpath(@__DIR__, "src", "genomic.jl"))


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Path to merged peaks bed"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Prefix of output file"
            arg_type = String
            required = true
        "bam"
            help = "Path to bam file"
            required = true
            action = :store_arg
            nargs = '+'
    end

    return parse_args(s)
end

args = parse_commandline()


function filter(record)::Bool
    #=
    Filter low quality reads
    =#

    if !BAM.ismapped(record)
        return false
    end

    auxdata = Dict(BAM.auxdata(record))
    if haskey(auxdata, "NH") && auxdata["NH"] > 1
        return false
    end

    if BAM.flag(record) & SAM.FLAG_QCFAIL > 0
        return false
    end

    if BAM.flag(record) & SAM.FLAG_PROPER_PAIR == 0
        return false
    end

    return true
end


function determine_strand(record)::String
    strand = "*"
    flag = BAM.flag(record)
    if flag & SAM.FLAG_READ1 > 0
        if flag & SAM.FLAG_REVERSE > 0
            strand = "-"
        else
            strand = "+"
        end
    elseif flag & SAM.FLAG_READ2 > 0
        if flag & SAM.FLAG_REVERSE > 0
            strand = "+"
        else
            strand = "-"
        end
    end
    return strand
end


function count_reads(bam::Vector, region::Genomic.BED)::Dict{String, Int}
    res = Dict() 

    for b = bam
        key = string(basename(Path(b)))
        key = replace(key, ".Aligned.sortedByCoord.out.bam" => "")
        reader = open(BAM.Reader, b, index=string(b, ".bai"))
        for record in eachoverlap(reader, region.Chrom, region.Start:region.End)
            if !filter(record)
                continue
            end

            # read strand
            strand = determine_strand(record)

            if strand != region.Strand
                continue
            end

            if BAM.leftposition(record) < region.End && BAM.rightposition(record) > region.Start

                res[key] = get(res, key, 0) + 1
            end
        end
    end
    return res
end


function load_bed(input_file::String)
    beds = []  # Dict{String, Vector}()
    open(input_file, "r") do r
        seekend(r)
        fileSize = position(r)

        seekstart(r)
        p = Progress(fileSize, 1)   # minimum update interval: 1 second
        while !eof(r)

            line = split(strip(readline(r)), "\t")
            
            if line[length(line)] == "NA"
                continue
            end

            try
                bic = parse(Float64, line[length(line)])
                if isinf(bic)
                    continue
                end
            catch e
                continue
            end

            site = split(line[1], r"[:-]")

            strand = site[length(site)]
            if strand != "+"
                strand = "-"
            end

            chrom, start_pos, end_pos = site[1], parse(Int, site[2]), parse(Int, site[3])
            alpha = split(line[5], ",")

            if length(alpha) > 0
                try
                    for x = alpha
                        if x != ""
                            s = round(Int, parse(Float64, x))
                            
                            if strand == "+"
                                push!(beds, Genomic.BED(chrom, s, s + 1, line[1], string(length(beds) + 1), strand))
                            else
                                push!(beds, Genomic.BED(chrom, s - 1, s, line[1], string(length(beds) + 1), strand))
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

function load_utr(utr::String)::Vector
    res =  Vector()

    open(utr, "r") do r
        seekend(r)
        fileSize = position(r)

        seekstart(r)
        p = Progress(fileSize, 1)   # minimum update interval: 1 second
        while !eof(r)
            line = split(strip(readline(r)), "\t")

            if occursin("|", line[4])
                continue
            end
            b = Genomic.BED(
                line[1],  parse(Int, line[2]), parse(Int, line[3]),
                line[4], line[5],  line[6]
            )

            push!(res, b)                          

            update!(p, position(r))
        end
        close(r)
    end

    return res
end


function main(input_file::String, bam::Vector, output::String)
    beds = load_utr(input_file)

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

    colnames = Vector()

    for b = bam
        key = string(basename(Path(b)))
        key = replace(key, ".Aligned.sortedByCoord.out.bam" => "")
        push!(colnames, key)
    end

    open(output, "w+") do w
        stream = nothing

        if endswith(output, ".gz")
            stream = ZlibDeflateOutputStream(w)
        end

        write(isnothing(stream) ? w : stream, string(",", join(colnames, ","), "\n"))

        p = Progress(length(beds), 1, "Computing...")
        Threads.@threads for i = eachindex(beds)
            b = beds[i]
            res = count_reads(bam,  b)

            row = [Genomic.get_bed_short(b)]
            for col = colnames
                push!(row, string(get(res, col, 0)))
            end
            write(isnothing(stream) ? w : stream, string(join(row, ","), "\n"))
            next!(p)

            b = nothing
            row = nothing
            res = nothing
            
            GC.safepoint()

            if i % 100 ==  0
                GC.gc()
            end
        end

        if !isnothing(stream)
            close(stream)
        end

        close(w)
    end
end


main(args["input"], args["bam"], args["output"])
