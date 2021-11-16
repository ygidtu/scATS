#!/usr/bin/env julia
#=
Created by Zhang Yiming at 2020.10.15

I try to convert apamix from python to julia
=#
using ArgParse
using FilePathsBase
using Memento
using ProgressMeter


include(joinpath(@__DIR__, "src", "genomic.jl"))


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--gtf", "-g"
            help = "Path to gtf file"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Path to output file"
            arg_type = String
            required = true
        "--expand", "-e"
            help = "The distance that utr to expand."
            arg_type = Int
            default = 500
    end

    return parse_args(s)
end


function call(cmd::Cmd, output::String)
    out = read(cmd, String)
    write(output, out)  
end


function main()
    parsed_args = parse_commandline()

    expand = parsed_args["expand"]
    gtf = parsed_args["gtf"]
    output =  parsed_args["output"]

    logger = Memento.config!("info"; fmt="[{date} | {level} | {name}]: {msg}")
    debug(logger, string(parsed_args))
    gtf = absolute(Path(gtf))
    if !exists(gtf)
        error(logger, string(gtf, " not exists"))
        exit(1)
    elseif !isfile(gtf)
        error(logger, string(gtf, " is not file"))
        exit(1)
    end

    output = absolute(Path(output))
    out_dir = parent(output)
    if !exists(out_dir)
        mkdir(out_dir)    
    end

    info(logger, string("reading from ", gtf))
    data = Genomic.load_GTF(string(gtf))
    
    # iterover transcripts
    utr_lst = Dict{String, Vector{Genomic.BED}}()

    info(logger, string("iter over all transcripts: ", length(data["transcript"])))
    @showprogress for (tid, transcript) in data["transcript"]
        if !haskey(data["exon"], tid)
            continue
        end
        exons = data["exon"][tid]
        exons = sort(exons)

        if transcript.Strand == "-"
            en = exons[length(exons)]
        else
            en = exons[1]
        end

        site = en.Strand == "+" ? en.Start : en.End

        # collect the exon closest to utr region
        bed = Genomic.BED(
            Chrom = en.Chrom, Start = max(1, en.Start - expand), End = en.End + expand,
            # Start = site - expand,
            # End = site + expand,
            Name = join([
                get(en.Attributes, "gene_biotype", "unknown"), 
                en.GeneID, en.TranscriptID, 
                get(en.Attributes, "gene_name", "unknown")
            ], ";"),
            Score = ".", Strand = en.Strand
        )

        key = string(bed.Chrom, "#", bed.Strand)

        temp = get(utr_lst, key, Vector())
        push!(temp, bed)
        utr_lst[key] = temp
    end

    info(logger, "write utr")
    
    res = Vector()
    for lst = values(utr_lst)
        append!(res, Genomic.merge(lst))
    end

    open(output, "w+") do w
        for i = sort(res)
            write(w, string(Genomic.get_bed(i), '\n'))
        end

        close(w)
    end
end

main()


#=
src = ""
o = ""
with open(o, "w+") as w:
    with open(src) as r:
        for line in r:
            line = line.split()
            site = line[1].split(":")
            pos = site[1].split("-")

            w.write(f"{site[0]}\t{pos[0]}\t{pos[1]}\t{line[0]}\t{line[1]}\t{site[2]}\n")
=#
