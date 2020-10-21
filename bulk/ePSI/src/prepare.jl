#!/usr/bin/env julia
#=
Created by Zhang Yiming at 2020.10.15

I try to convert apamix from python to julia
=#
using ArgParse
using FilePathsBase
using Formatting
using Memento


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
        "--aggregate", "-a"
            help = "Indicates whether two or more genes sharing an exon should be merged into an 'aggregate gene'. If 'no', the exons that can not be assiged to a single gene are ignored."
            arg_type = Bool
            default = false
    end

    return parse_args(s)
end


function call(cmd::Cmd, output::String)
    out = read(cmd, String)
    write(output, out)  
end

function main(gtf::String, output::String, aggregate::Bool)
    logger = Memento.config!("info"; fmt="[{date} | {level} | {name}]: {msg}")
    
    gtf = absolute(Path(gtf))
    if !exists(gtf)
        error(logger, string(gtf, " not exists"))
        exit(1)
    end

    output = absolute(Path(output))
    out_dir = parent(output)
    if !exists(out_dir)
        mkdir(out_dir)    
    end

    curr_dir = @__DIR__
    cmd = `python $curr_dir/prepare.py -i $gtf -o stdout`
    if aggregate
        cmd = `python $curr_dir/prepare.py -i $gtf -o stdout -r`
    end
    
    fe = "{}\t{}:{}\n"
    open(output, "w+") do w
        open(cmd) do io
            while !eof(io)
                lines = split(strip(replace(readline(io), r"[\";]" => "")), r"\s")

                if lines[3] != "exonic_part"
                    continue
                end

                write(w, Formatting.format(fe, join(lines[1:8], "\t"), lines[14], lines[12]))
            end
            close(io)
        end
        close(w)
    end
end

args = parse_commandline()
main(args["gtf"], args["output"], args["aggregate"])
