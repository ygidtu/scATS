using ArgParse
using FilePathsBase
using Memento
using ProgressMeter

#=
using StatsBase
=#
include(joinpath(@__DIR__, "../modeling/src", "genomic.jl"))

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Path to input file"
            arg_type = String
            required = true
        "--cage", "-c"
            help = "Path to cage bed"
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

logger = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")


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


function main(input_file::String, cage::String, output::String)
    info(logger, "load ATS")
    beds = load_bed(input_file)

    temp = string(output, ".ats.bed")
    open(temp, "w+") do w
        for b = beds
            write(w, string(b, "\n"))
        end
    end

    temp_cage = string(output, ".temp_cage")
    open(temp_cage, "w+") do w
        open(cage, "r") do r
            while !eof(r)
                line = readline(r)

                write(w, string(replace(line, r"^chr" => ""), "\n"))
            end
            close(r)
        end
        close(w)
    end

    info(logger, "dist")
    closest_features = joinpath(@__DIR__, "bin/closest-features")
    run(pipeline(`$closest_features --closest --dist $temp $temp_cage`, stdout=output))

    if exists(Path(temp_cage))
        rm(Path(temp_cage))
    end
end


main(args["input"], args["cage"], args["output"])

