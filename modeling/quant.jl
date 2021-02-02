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
    using BioAlignments
    using FilePathsBase
    using GenomicFeatures
    using Memento
    using ProgressMeter
    using XAM
    
    #=
    using StatsBase
    =#
    include(joinpath(@__DIR__, "src", "genomic.jl"))

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

    function count_reads(bam::String, region::Genomic.BED, barcodes::Vector{String})::Dict{String, Int}
        res = Dict() 

        reader = open(BAM.Reader, bam, index=string(bam, ".bai"))
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
                auxdata = Dict(BAM.auxdata(record))
                cb = auxdata["CB"]
                ub = auxdata["UB"]

                if !haskey(res, cb)
                    res[cb] = Dict()
                end
    
                res[cb][ub] = get(res[cb], ub, 0) + 1
            end
        end

        # key = collect(keys(res))
        # for k = key[1:5]
        #     println(string("key: ", k, " ", res[k]))
        # end

        return Dict(x => length(y) for (x, y) = res)
    end
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

            # if length(beds) >= 5
            #     break
            # end

            update!(p, position(r))
        end
        close(r)
    end

    return beds
end


function main(input_file::String, cellranger::String, output::String)
    beds = load_bed(input_file)

    # read barcodes
    barcodes = Vector{String}()
    barcode = joinpath(cellranger, "filtered_feature_bc_matrix/barcodes.tsv.gz")
    for line in eachline(open(barcode) |> ZlibInflateInputStream)
        push!(barcodes, line)
    end

    bam = joinpath(cellranger, "possorted_genome_bam.bam")

    res = @showprogress 1 "Counting..." pmap(beds) do b
        res = count_reads(bam,  b, barcodes)

        row = [Genomic.get_bed_short(b)]
        for barcode = barcodes
            push!(row, string(get(res, barcode, get(res, replace(barcode, "-1"=> ""), 0))))
        end

        return row
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
    output = string(output)
    open(output, "w+") do w

        stream = nothing

        if endswith(output, ".gz")
            stream = ZlibDeflateOutputStream(w)
        end

        write(isnothing(stream) ? w : stream, string(",", join(barcodes, ","), "\n"))

        for row = res
            write(isnothing(stream) ? w : stream, string(join(row, ","), "\n"))
        end

        if !isnothing(stream)
            close(stream)
        end

        close(w)
    end
end


main(args["input"], args["cellranger"], args["output"])
