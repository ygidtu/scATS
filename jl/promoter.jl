#!/usr/bin/env julia
#=
Created by Zhang Yiming at 2020.12.22

Using R proActiv pacakges to prepare the promoter region
=#
using ArgParse
using FilePathsBase
using Memento
using ProgressMeter



function main(gtf::String, output::String; species::String = "Homo_sapiens")
    rcopy(R"""
    library(proActiv)
    promoterAnnotation <- preparePromoterAnnotation(
        file=$gtf,
        species = species
    )

    temp = as.data.frame(promoterAnnotation@promoterCoordinates)
    temp1 = as.data.frame(promoterAnnotation@promoterIdMapping)
    temp = merge(temp, temp1, by = "promoterId")

    write.csv(temp, $output, quote = F, row.names = F)
    """)

end



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--gtf", "-g"
            help = "Path to gtf file"
            arg_type = String
            required = true
        "--species", "-s"
            help = "Species of gtf file"
            arg_type = String
            default="Homo_sapiens"
        "--output", "-o"
            help = "Output file path"
            arg_type = String
    end

    return parse_args(s)
end


parsed_args = parse_commandline()

main(parsed_args["gtf"], parsed_args["output"], species=parsed_args["species"])