#!/usr/bin/env julia
#=
Created by Zhang Yiming at 2020.10.27
=#
using ArgParse
using Distributed


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input-bam", "-i"
            help = "Path to gtf file"
            arg_type = String
            required = true
        "--process", "-p"
            help = "How many CPU to  use"
            arg_type = Int64
            default = 1
        "--output", "-o"
            help = "Path to output file"
            arg_type = String
            required = true
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
    using Formatting
    using GenomicFeatures
    using ProgressMeter

    struct BEDData
        chrom::String
        start_pos::Int64
        end_pos::Int64
        strand::String
        score::String
        name::String 
    end

    Base.show(io::IO, self::BEDData) = print(
        io,
        Formatting.format(
            FormatExpr("{}:{}-{}:{}\t{}\t{}"),
            self.chrom, self.start_pos, self.end_pos, self.strand,
            self.name, self.score
        )
    )
    
    function new_bed(line::String)::BEDData
        lines = split(strip(line), "\t")
    
        return BEDData(
            lines[1], parse(Int64, lines[2]), parse(Int64, lines[3]),
            lines[6], lines[4], lines[5]
        )
    end

    function filter(record)::Bool
        #=
        Filter low quality reads
        =#

        auxdata = Dict(BAM.auxdata(record))
        if haskey(auxdata, "NH") && auxdata["NH"] > 1
            return false
        end

        if BAM.flag(record) & SAM.FLAG_DUP > 0 || BAM.flag(record) & SAM.FLAG_QCFAIL > 0
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

    function process_bam(path::String, ref::String, size::Int64)::Vector{String}
        res = Vector{String}()

        start_sites = Dict{String, Dict{Int, Int}}()

        # println(string("reader: ", reader, "; typeof: ", typeof(reader)))
        # println(string("ref: ", ref, "; typeof: ", typeof(ref)))
        # println(string("size: ", size, "; size: ", typeof(size)))
        bai = string(path, ".bai")
        open(BAM.Reader, path, index=bai) do reader
            for record in eachoverlap(reader, Interval(ref, 1, size))
                if !filter(record)
                    continue
                end

                strand = determine_strand(record)
                if strand == "*"
                    continue
                end
                start_pos = BAM.position(record)

                if strand == "-"
                    start_pos = BAM.rightposition(record)
                end

                cb = BAM.auxdata(record)["CB"]
                ub = BAM.auxdata(record)["UB"]
                cb = string(cb, "\t", ub)
                temp_sites = get(start_sites, cb, Dict())
                temp_sites[start_pos] = get(temp_sites, start_pos, 0) + 1
                start_sites[cb] = temp_sites
            end
            close(reader)
        end

        for (cb, sites) in start_sites
            temp = [Format.format(FormatExpr("{}:{}"), i, j) for (i, j) in sites]
            push!(
                res, 
                Formatting.format(
                    FormatExpr("{}\t{}"),
                    cb, join(sort(temp), ";")
                )
            )
        end

        return res
    end
end


function references(reader)::Dict{String, Int64}
    h = header(reader)
    res = Dict()

    for sq in findall(h, "SQ")
        vals = values(sq)
        res[vals[1]] = parse(Int64, vals[2])
    end
    return res
end


function main(input_bam::String, output::String)
    out_dir = parent(absolute(Path(output)))
    if !exists(out_dir)
        mkdir(out_dir, recursive=true)
    end
    # Load genomic features from a BED file.
    bai = string(input_bam, ".bai")
    if !exists(Path(bai))
        run(`samtools index $bai`)
    end

    # 1:1-248956422
    reader = open(BAM.Reader, input_bam, index=bai)
    refs = references(reader)
    close(reader)

    res = @showprogress pmap(collect(keys(refs))) do ref
        try
            return process_bam(input_bam, ref, refs[ref])
        catch e
            println(string(ref, ": ", e))
            return []
        end
    end

    open(output, "w+") do w
        for r in res
            write(w, join(r, "\n"))
            write(w, "\n")
        end
        close(w)
    end
end

main(args["input-bam"], args["output"])


