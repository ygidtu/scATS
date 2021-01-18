#!/usr/bin/env julia
#=
Created at 2021.01.18

This script is used to quantify the expression counts of each peaks
=#

using ArgParse
using Distributed
using Libz


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Path to merged peaks bed"
            arg_type = String
            required = true
        "--cellranger", "-c"
            help = "Path to cellranger outs directory"
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
    
    #=
    using StatsBase
    =#
    include(joinpath(@__DIR__, "src", "extract.jl"))

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

        return haskey(auxdata, "CB") && haskey(auxdata, "UB")
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

    function count_reads(bam::String, region::Extract.BEDRecord, barcodes::Vector{String})::Dict{String, Int}
        res = Dict(x => Dict() for x = barcodes)

        reader = open(BAM.Reader, path, index=string(bam, ".bai"))

        for record in eachoverlap(reader, region.chrom, region.start_pos:region.end_pos)
            if !filter(record)
                continue
            end

            # read strand
            strand = determine_strand(record)

            if strand != region.strand
                continue
            end
            auxdata = Dict(BAM.auxdata(record))
            cb = auxdata["CB"]
            ub = auxdata["UB"]

            res[cb][ub] = 1
        end

        return Dict(x => length(y) for (x, y) = res)
    end
end


function main(input_file::String, cellranger::String, output::String)
    beds = []  # Dict{String, Vector}()
    open(input_file, "r") do r
        seekend(r)
        fileSize = position(r)

        if verbose
            Extract.setLevel("debug")
        end

        seekstart(r)
        p = Progress(fileSize, 1)   # minimum update interval: 1 second
        while !eof(r)
            temp_bed = Extract.new_bed(readline(r)) #, chromosomes=chromosomes)

            push!(beds, temp_bed)

            if length(beds) > 500
                break
            end

            update!(p, position(r))
        end
        close(r)
    end

    # read barcodes
    barcodes = Vector{String}()
    barcode = joinpath(cellranger, "filtered_feature_bc_matrix/barcodes.tsv.gz")
    for line in eachline(open(barcode) |> ZlibInflateInputStream)
        push!(barcodes, replace(line, "-1"=>""))
    end

    bam = joinpath(cellranger, "possorted_genome_bam.bam")

    open(output, "w+") do w

        stream = nothing

        if endswith(output, ".gz")
            stream = ZlibDeflateOutputStream(w)
        end

        write(w, string(",", ",".join(barcodes), "\n"))

        @showprogress 1 "Counting..." pmap(beds) do b
            res = count_reads(bam,  b, barcodes)

            row = [Extract.get_bed_short(b)]
            for barcode = barcodes
                push!(row, string(res[barcode]))
            end

            write(isnothing(stream) ? w : stream, string(",".join(row), "\n"))
        end

        if !isnothing(w)
            close(stream)
        end

        close(w)
    end
end


main(args["input"], args["cellranger"], args["output"])
