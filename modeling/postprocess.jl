#!/usr/bin/env julia
#=
Created at 2021.01.18

This script is used to quantify the expression counts of each peaks
=#

using ArgParse
using FilePathsBase
using Libz
using Memento
using ProgressMeter


#=
using StatsBase
=#
include(joinpath(@__DIR__, "src", "genomic.jl"))


function load_bed(input_file::String)
    beds = Dict{String, String}()

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
                                b = Genomic.BED(chrom, s, s + 1, line[1], string(length(beds) + 1), strand)

                                beds[Genomic.get_bed_short(b)] = line[1]
                            else
                                b = Genomic.BED(chrom, s - 1, s, line[1], string(length(beds) + 1), strand)

                                beds[Genomic.get_bed_short(b)] = line[1]
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


function load_utr(utr::String)::Dict{String, String}
    res =  Dict{String, String}()

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

            n = split(line[4], ";")
            n = n[length(n)]

            res[Genomic.get_bed_short(b)] = n                            

            update!(p, position(r))
        end
        close(r)
    end

    return res
end


function main(input_file::String, utr::String, output::String)

    beds = load_bed(input_file)
    utr = load_utr(utr)

    open(output, "w+") do w
        for (key, val) = beds
            if haskey(utr, val)
                write(w, string(val, "\t", key, "\t", utr[val], "\n"))
            end
        end
    end
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Path to merged peaks bed"
            arg_type = String
            required = true
        "--utr", "-u"
            help = "Path to utr bed"
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


main(args["input"], args["utr"], args["output"])
