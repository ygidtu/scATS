#!/usr/bin/env julia
#=
Created by Zhang Yiming at 2020.10.15

I try to convert apamix from python to julia
=#
using ArgParse
using FilePathsBase
using Memento
using ProgressMeter

include("genomic.jl")


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--gtf", "-g"
            help = "Path to gtf file"
            arg_type = String
            required = true
        "--prefix", "-p"
            help = "Prefix of output file"
            arg_type = String
            required = true
        "--distance", "-d"
            help = "The distance between utr to merge."
            arg_type = Int
            default = 0
    end

    return parse_args(s)
end


function call(cmd::Cmd, output::String)
    out = read(cmd, String)
    write(output, out)  
end


function main()
    parsed_args = parse_commandline()

    distance = parsed_args["distance"]

    logger = Memento.config!("info"; fmt="[{date} | {level} | {name}]: {msg}")
    debug(logger, string(parsed_args))
    gtf = absolute(Path(parsed_args["gtf"]))
    if !exists(gtf)
        error(logger, string(gtf, " not exists"))
        exit(1)
    end

    output = absolute(Path(parsed_args["prefix"]))
    out_dir = parent(output)
    if !exists(out_dir)
        mkdir(out_dir)    
    end

    data = Dict(
        "transcript"=>Dict{String, Genomic.GTF}(), 
        "exon"=>Dict{String, Vector}()
    )

    info(logger, string("reading from ", gtf))
    open(gtf) do io
        seekend(io)
        fileSize = position(io)
        seekstart(io)

        p = Progress(fileSize, 1)   # minimum update interval: 1 second
        while !eof(io)
            line = readline(io)
            if startswith(line, "#")
                continue
            end

            gtf = Genomic.create_GTF(line)
            if gtf.Type == "transcript"
                # if there is not gene_biotype or gene_biotype is not specific types passed
                if !haskey(gtf.Attributes, "gene_biotype") || !(gtf.Attributes["gene_biotype"] in ["antisense", "lincRNA", "protein_coding"])
                    continue
                end

                data[gtf.Type][gtf.ID] = gtf
            elseif gtf.Type == "exon"
                if !haskey(data["transcript"], gtf.TranscriptID)
                    continue
                end

                if !haskey(data[gtf.Type], gtf.TranscriptID)
                   data[gtf.Type][gtf.TranscriptID] = Vector()
                end
                push!(data[gtf.Type][gtf.TranscriptID], gtf)
            end
            update!(p, position(io))
        end
    end
    
    # iterover transcripts
    mt_lst = Vector{Genomic.BED}()
    utr_lst = Vector{Genomic.BED}()
    intron_lst = Vector{Genomic.BED}()
    exon_lst = Vector{Genomic.BED}()
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

        # collect the exon closest to utr region
        bed = Genomic.BED(
            en.Chrom, en.Start, en.End,
            join([
                get(en.Attributes, "gene_biotype", "unknown"), 
                en.GeneID, en.TranscriptID, 
                get(en.Attributes, "gene_name", "unknown")
            ], ";"),
            ".", en.Strand
        )

        if uppercase(en.Chrom) == "MT"
            push!(mt_lst, bed)
        else
            push!(utr_lst, bed)
        end

        # collect the intron
        for i = 2:2:length(exons)
            bed = Genomic.BED(
                exons[i].Chrom,
                exons[i-1].End + 1,
                exons[i].Start - 1,
                join([
                    exons[i].Chrom, exons[i].GeneID, exons[i].TranscriptID,
                    get(exons[i].Attributes, "gene_name", "unknown")
                ], ";"),
                ".", exons[i].Strand
            )
            push!(intron_lst, bed)

            if i == 2
                bed = Genomic.BED(
                    exons[i-1].Chrom,
                    exons[i-1].Start,
                    exons[i-1].End,
                    join([
                        exons[i-1].Chrom, exons[i-1].GeneID, exons[i-1].TranscriptID,
                        get(exons[i-1].Attributes, "gene_name", "unknown")
                    ], ";"),
                    ".", exons[i-1].Strand
                )
                push!(exon_lst, bed)
            end
            
            bed = Genomic.BED(
                exons[i].Chrom,
                exons[i].Start,
                exons[i].End,
                join([
                    exons[i].Chrom, exons[i].GeneID, exons[i].TranscriptID,
                    get(exons[i].Attributes, "gene_name", "unknown")
                ], ";"),
                ".", exons[i].Strand
            )
            push!(exon_lst, bed)
        end
    end

    info(logger, "write utr")
    utr = string(output, "_utr.bed")
    temp_utr = string(utr, ".tmp")
    open(temp_utr, "w+") do w
        append!(mt_lst, utr_lst)
        for i = sort(mt_lst)
            write(w, string(i, '\n'))
        end
        close(w)
    end

    call(`bedtools merge -i $temp_utr`, utr)
    rm(temp_utr)

    exon_bed = string(output, "_exon.bed")
    open(exon_bed, "w+") do w
        for i = sort(exon_lst)
            write(w, string(i, '\n'))
        end
    end

    info(logger, "write intron")
    intron_bed = string(output,  "_intron.bed")
    temp_intron = string(intron_bed, ".tmp")
    open(temp_intron, "w+") do w
        intron_lst = Genomic.merge(intron_lst)
        for i = sort(intron_lst)
            write(w, string(i, '\n'))
        end
        close(w)
    end
    call(`bedtools subtract -a $temp_intron -b $exon_bed`, intron_bed)
    rm(temp_intron)
end

main()

