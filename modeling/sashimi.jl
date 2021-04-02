
using ArgParse
using Distributed
# using RCall


function parse_commandline()
    s = ArgParseSettings("pipeline to make venn diagram using R VennDiagram package")

    @add_arg_table! s begin
        "--input", "-i"
            help = "Path to QTL output directory"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Path to output file"
            arg_type = String
            required = true
        "--bam", "-b"
            help = "Path to STAR output directory"
            arg_type = String
            required = true
        "--reference", "-r"
            help = "Path to reference gtf"
            arg_type = String
            required = true
        "--process", "-p"
            help = "How many processes to use"
            arg_type = Int64
            default = 10
        "--reverse"
            help = "plot the sites with biggest BIC"
            default = false
            action = :store_true
    end

    return parse_args(s)
end


args = parse_commandline()

if get(args, "process", 1) > 1
    addprocs(get(args, "process", 1))
end


@everywhere begin
    using FilePathsBase
    using ProgressMeter
    
    mutable struct Site
        Chrom::String
        Start::Int64
        End::Int64
        Strand::String
        Counts::Int64
        BIC::Float64
        Site::String
        RealST::String
        RealEn::String
    end

    Base.show(io::IO, self::Site) = print(io, string(self.Chrom, ":", self.Start, "-", self.End, ":", self.Strand, "_", round(self.BIC, digits = 2), "_", self.Counts))

    struct Bam
        CellType::String
        Color::String
        ID::String
        Path::String
    end

    Base.show(io::IO, self::Bam) = print(io, join([self.Path, self.CellType, self.Color, self.ID], "\t"))


    function decode_attrs(str::AbstractString)::Dict{String, String}
        str = replace(str, "\""=>"")
        lines = split(str, ";")
    
        res = Dict()
        for i = lines
            j = split(strip(i), " ")
            if length(j) < 2
                continue
            end
            res[strip(j[1])] = strip(j[2])
        end
        return res
    end
    

    function get_gene_region(ref::String, site::Site)::Dict{String, Site}
        res = Dict()

        region = string(site.Chrom, ":", site.Start, "-", site.End)
        open(`tabix $ref $region`) do r
            for line in eachline(r)
                line = split(strip(line), "\t")
                if line[3] == "gene"
                    attr = decode_attrs(line[9])
                    bic = site.BIC
                    if bic < 0
                        bic *= -1
                    end
                    res[get(attr, "gene_id", "")] = Site(line[1], parse(Int64, line[4]), parse(Int64, line[5]), line[7], site.Counts, bic)
                end
            end
            close(r)
        end
        return res
    end


    function plot(bam::String, site::Site, ref::String, output::String, script::String, process::Int)

        # for (gene, region) in get_gene_region(replace(ref, ".gtf"=>".sorted.gtf.gz"), site)
    
        # dist = round(Int, (site.End - site.Start) / 2)
        # dist = 0
        dist = round(Int, site.End - site.Start)

        # lines = join([string(site.Start), string(site.End)], ",")
        lines = site.Site
        event = string(site.Chrom, ":", max(site.Start - dist, 1), ":", site.End + dist) # , ":", site.Strand
        focus = string(site.Start, "-", site.End)
        dots = site.RealST
        dots2 = site.RealEn

        o = joinpath(output, string(site, ".pdf"))

        # if string(site) != "1:19596110-19597432:-_76929.83_2"
        #     return
        # end

        # /mnt/raid61/Personal_data/zhangyiming/software/pyenv/versions/3.7.6/bin/
        run(`sashimiplot junc --gtf $ref --bam $bam --sj 100 --junc $event --ie 1,1  --ps RF --ssm R1 --fileout $o --trackline $lines --focus $focus`)
        # run(`python $script normal -b $bam -e $event --indicator-lines $lines -g $ref -o $o --log 0 -t 15 -p $process`)
        
    end
end


function load(path::String; max_size::Int64=100, reverse::Bool = false, counts::Int=20)::Vector{Site}

    res = Vector{Site}()
    r = open(path)

    seekend(r)
    fileSize = position(r)

    seekstart(r)
    p = Progress(fileSize, 1)   # minimum update interval: 1 second
    while !eof(r)
        line = split(strip(readline(r)), "\t")
        
        # reads_starts = parse.(Int, split(line[2], ","))
        # reads_ends = parse.(Int, split(line[3], ","))
        # if length(reads_starts) < counts
        #     continue
        # end

        if line[length(line)] == "NA"
            continue
        end

        try
            bic = parse(Float64, line[length(line)])

            if isinf(bic)
                continue
            end
        catch e
            # println(e)
            # println(line[length(line)])
            continue
        end

        if reverse
            bic *= -1
        end

        site = split(line[1], r"[:-]")

        strand = site[length(site)]
        if strand != "+"
            strand = "-"
            # continue
        end
        # println(site)
        start_pos, end_pos = parse(Int, site[2]), parse(Int, site[3])
        lab_sites = [site[4] == "+" ? end_pos : start_pos]
        # println(line[5])
        alpha = split(line[3], ",")
        ws = parse.(Float64, split(line[2], ","))
        # beta = split(line[6], ",")

        alpha = [y for (x, y) = zip(ws, alpha) if x > 0.4]

        if length(alpha) < 2
            continue
        end

        try
            for x = alpha
                if x != ""
                    s = site[4] != "+" ? start_pos + round(Int, parse(Float64, x)) : end_pos - round(Int, parse(Float64, x))
                    push!(lab_sites, s)

                    # push!(lab_sites, round(Int, parse(Float64, x)))

                end
            end
        catch e
            println(e)
            continue
        end
        
        if length(lab_sites) < 1
            continue
        end

        site = Site(
            site[1], 
            start_pos, end_pos,
            # reads_starts[1], 
            # reads_ends[length(reads_ends)],
            strand,
            length(split(line[2], ",")), 
            bic,
            join(map(string, lab_sites), ","),
            line[4],
            line[5]
        )
        push!(res, site)
        # if length(res) < 3
        #     push!(res, site)
        #     res = sort!(res, by=x -> x.BIC)
        # elseif res[length(res)].BIC < bic
        #     push!(res, site)
        # else
        #     low, high = 1, length(res)
        #     # println("x=", x, "; i=", i)
        #     while high > low + 1
        #         mid = floor(Int, (low + high) / 2)
        #         # println(string("low: ", low, "; mid=", mid, "; high=", high))
        #         if res[mid].BIC < bic
        #             low = mid
        #         elseif res[mid].BIC > bic
        #             high = mid
        #         else
        #             # low = mid
        #             break
        #         end  
        #     end
            
        #     insert!(res, low + 1, site)

        #     if length(res) > max_size
        #         res = res[1:max_size]
        #     end
        # end

        update!(p, position(r))
    end
    close(r)

    return res
end


function main(input::String; bam::String, ref::String, script::String, output::String, reverse::Bool=false, process::Int=1)
    if !exists(Path(output))
        mkdir(Path(output), recursive=true)
    end
    data = load(input, max_size=100)
    @showprogress pmap(data) do g
        plot(bam, g, ref, output, script, process)
    end
end

# input = "/mnt/raid64/Chen_Cell_2016/alignment/STAR/"
# rds = "rds/PSI_mQTL.rds"
# ref =  "ref/gencode.v30lift37.annotation.exonic_part.sorted.gff"
# "/mnt/raid61/Personal_data/zhangyiming/ePSI/rds/QTL"
# "/mnt/raid61/Personal_data/zhangyiming/ePSI/rds/QTL/venn.pdf"

main(
    args["input"], 
    bam=args["bam"], 
    ref=args["reference"],
    script="/mnt/raid61/Personal_data/zhangyiming/code/pysashimi/main.py",
    output=args["output"],
    reverse=args["reverse"],
    process=args["process"]
)

