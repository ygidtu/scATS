#!/usr/bin/env julia
#=
Created by Zhang Yiming at 2020.11.27

This script is used to merge all infered start sites
=#
using ArgParse
using FilePathsBase
using Parameters
using ProgressMeter


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--output", "-o"
            help = "Path to output."
            arg_type = String
            required = true
        "--expand", "-e"
            help = "How many bp to expand"
            arg_type = Int
            default = 100
        "ats"
            help = "Path to atsmix output file"
            required = true
            action = :store_arg
            nargs = '+'
    end

    return parse_args(s)
end


@with_kw struct Site
    Chrom::String
    Start::Int64
    Strand::String
    ID::String
end

Base.show(io::IO, self::Site) = print(io, Formatting.format(FormatExpr("{}:{}:{}"), self.chrom, self.Start, self.Strand))


function get_site(site::Site, expand::Int=100)::String
    half = round(Int, expand / 2)
    return Formatting.format(
        FormatExpr("{}\t{}\t{}\t.\t.\t{}"), 
        self.chrom, max(1, site.Start - half), 
        site.Start + half, self.Strand
    )
end


function load(paths::Vector{String}, distance::Int64=100)::Dict{String, Vector{Site}}
    # res contains the site id -> list of site
    res = Dict{String, Vector{Site}}()

    @showprogress 1 "Reading..." for path = paths
        r = open(path)

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

            start_pos, end_pos = parse(Int, site[2]), parse(Int, site[3])
            alpha = split(line[5], ",")

            if length(alpha) > 0
                for x = alpha
                    if x != ""
                        s = site[4] != "+" ? start_pos + round(Int, parse(Float64, x)) : end_pos - round(Int, parse(Float64, x))

                        id = string(site[1], ":", round(s / distance))
                        site = Site(site[1], s, strand, id)

                        temp = get(res, id, Vector())
                        append!(temp, site)
                        res[id] = temp
                    end
                end
            end
            update!(p, position(r))
        end
        close(r)
    end

    return res
end


function write_to(data::Dict{String, Vector{Site}}, output::String; expand::Int=100)

    open(output, "w+") do w

        for res = collect(values(data))
            sites = Dict{String, Int}()
            for r = res
                site = get_site(r, expand)
                sites[site] = get(sites, site, 0) + 1
            end
            most_expr_site = argmax(collect(values(sites)))

            if most_expr_site > 1
                site = collect(keys(sites))[filter(x -> x == most_expr_site, collect(values(sites)))]

                if length(site) >  0
                    write(w, string(site, "\n"))
                end
            else
                pos = sort([x.Start for x = res])
                write(w, Formatting.format(
                    FormatExpr("{}\t{}\t{}\t.\t.\t{}\n"), 
                    res[1].Chrom, pos[1], 
                    pos[length(pos)], res[1].Strand
                ))
            end
        end
            
        close(w)
    end
end


function main(
    ats::Vector{String}, 
    output::AbstractString; 
    expand::Int=100
)
    write_to(load(ats, expand), output; expand=expand)
end


parsed_args = parse_commandline()
println(parsed_args)
main(
    parsed_args["ats"], 
    parsed_args["output"],
    expand=parsed_args["expand"]
)
