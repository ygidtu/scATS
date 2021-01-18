#!/usr/bin/env julia

using ArgParse
using BioGenerics
using Distributed
using XAM


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Path to utr bed"
            arg_type = String
            required = true
        "--bam", "-b"
            help = "Path to bam file"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Prefix of output file"
            arg_type = String
            required = true
        "--distance", "-d"
            help = "The minimum distance of read in utr."
            arg_type = Int
            default = 1500
        "--using-R"
            help="whether to use R ksmooth"
            action = :store_true
        "--process", "-p"
            help = "How many processes to use"
            arg_type = Int64
            default = 1
        "--n-max-ats"
            help = "the maximum of ats inside a utr"
            arg_type = Int64
            default = 5
        "--n-min-ats"
            help = "the minimum of ats inside a utr"
            arg_type = Int64
            default = 1
        "--mu-f"
            help = "mu f"
            arg_type = Int64
            default = 300
        "--sigma-f"
            help = "sigma-f"
            arg_type = Int64
            default = 50
        "--min-ws"
            help = "min ws"
            arg_type = Float64
            default = 0.01
        "--max-beta"
            help = "maximum beta"
            arg_type = Float64
            default = 50.0
        "--min-reads"
            help = "minimum reads to construct ATS"
            arg_type = Int
            default = 5
        "--seed"
            help = "seed"
            arg_type = Int64
            default = 42
        "--single-end"
            help="whether this is sinle-end sequencing"
            action = :store_true
        "--fixed-inference"
            help="inference with fixed parameters"
            action = :store_true
        "--verbose"
            help="whether to display additional message"
            action = :store_true
    end

    return parse_args(s)
end

args = parse_commandline()
# println(args)
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
end


function main(
    input_file::String, bam::String, output::String, distance::Int;
    mu_f::Int64=300, min_ws::Float64=0.01,
    max_beta::Float64=50.0, sigma_f::Int64=50,
    using_R::Bool=true, verbose::Bool=false,
    process::Int=1, n_max_ats::Int=5, n_min_ats::Int=1,
    fixed_inference::Bool=false, single_end::Bool=false,
    min_reads::Int=0
)
    logger = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")
    if verbose
        logger = Memento.config!("debug"; fmt="[{date} - {level} | {name}]: {msg} | {stacktrace}")
    end
    output = absolute(Path(output))
    out_dir = parent(output)
    try
        if !exists(out_dir)
            mkdir(out_dir, recursive=true)
        end
    catch e
        error(logger, Formatting.format(FormatExpr("Error while create {}: {}"), output, e))
        exit(1)
    end

    beds = []  # Dict{String, Vector}()
    r = open(input_file, "r")
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

    info(logger, string("The number of utrs: ", length(beds)))

    info(logger, "check index")
    bai = string(bam, ".bai")
    if !exists(Path(bai))
        run(`samtools index $bam`)
    end
    
    info(logger, "start running")
    error_log = string(output, "_error.log")
    if exists(Path(error_log))
        rm(error_log)
    end

    res = @showprogress "Reading... " pmap(beds) do b
        temp = Extract.get_record_from_bam(
            bam, b, 
            distance=distance, 
            min_reads=min_reads
        )

        temp_res = Vector()
        for data = temp
            append!(
                temp_res,
                Extract.run(
                    data; n_max_ats = n_max_ats, n_min_ats = n_min_ats,
                    mu_f = mu_f, sigma_f = sigma_f, min_ws = min_ws, 
                    max_beta = max_beta, fixed_inference_flag = fixed_inference,
                    single_end_mode = single_end, verbose = verbose, 
                    using_R = using_R, error_log = error_log, density=string("/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_R/density/", Extract.get_utr(data), ".pdf")
                )
            )
        end
        return temp_res
    end

    # # extract reads and drop duplicates
    # temp_data = Vector{Extract.ExtractedData}()
    # for t = res
    #     for r = t
    #         push!(temp_data, r)
    #     end
    # end

    # res = @showprogress "Computing... "  pmap(collect(values(temp_data))) do data
    #     return Extract.run(
    #         data; n_max_ats = n_max_ats, n_min_ats = n_min_ats,
    #         mu_f = mu_f, sigma_f = sigma_f, min_ws = min_ws, 
    #         max_beta = max_beta, fixed_inference_flag = fixed_inference,
    #         single_end_mode = single_end, verbose = verbose, 
    #         using_R = using_R, error_log = error_log, density=string("/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_R/density/", Extract.get_utr(data), ".pdf")
    #     )
    # end

    # res = Vector()
    # for data = collect(values(temp_data))
    #     push!(res, Extract.run(
    #         data; n_max_ats = n_max_ats, n_min_ats = n_min_ats,
    #         mu_f = mu_f, sigma_f = sigma_f, min_ws = min_ws, 
    #         max_beta = max_beta, fixed_inference_flag = fixed_inference,
    #         single_end_mode = single_end, verbose = verbose, 
    #         using_R = using_R, error_log = error_log,
    #     ))
    # end

    open(output, "w+") do w
        write(w, "utr\tst_arr\ten_arr\tws\talpha_arr\tbeta_arr\tlb_arr\tlabel\tbic\n")
        for t = res
            for r = t
                if r != ""
                    write(w, string(r, "\n"))
                end
            end
        end
        close(w)
    end
end


main(
    args["input"], args["bam"], args["output"], 
    args["distance"], mu_f = args["mu-f"], 
    min_ws = args["min-ws"], sigma_f = args["sigma-f"],
    max_beta=args["max-beta"], # using_R=args["using-R"], 
    verbose=args["verbose"], process=args["process"], 
    n_max_ats=args["n-max-ats"], n_min_ats=args["n-min-ats"],
    fixed_inference=args["fixed-inference"], 
    single_end=args["single-end"], min_reads=args["min-reads"]
)


function test(mu_f::Int64=350, min_ws::Float64=0.01,
    max_beta::Float64=50.0, sigma_f::Int64=70,
    using_R::Bool=true, verbose::Bool=false,
    n_max_ats::Int=5, n_min_ats::Int=2, distance::Int=1500,
    fixed_inference::Bool=false, single_end::Bool=false,
    min_reads::Int=0)

    bam = "/mnt/raid64/ATS/Personal/zhangyiming/bams/NHC2.bam"

    
    bed = Extract.new_bed("1\t1471469\t1472189\t.\t.\t+")
    # "1\t6205375\t6206201\t.\t.\t+"
    # 1\t2556266\t2556833
    
    println(bed)
    
    temp = Extract.get_record_from_bam(
        bam, bed, 
        distance=distance, 
        min_reads=min_reads
    )
    println(temp[1])
    res = Vector()
    for seed = 1:100
        # seed = 4
        res = Extract.run(
            temp[1]; n_max_ats = n_max_ats, n_min_ats = n_min_ats,
            mu_f = mu_f, sigma_f = sigma_f, min_ws = min_ws, 
            max_beta = max_beta, fixed_inference_flag = fixed_inference,
            single_end_mode = single_end, verbose = verbose, 
            using_R = using_R, error_log = nothing, seed=seed,
            debug = true, density="/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_R/test.pdf"
        )
        # println(string("temp res: ", res))
        if length(res) > 0 && !isnothing(res[1].alpha_arr) && length(res[1].alpha_arr)  > 0
            println(string("seed: ", seed))
            break
        end
    end
    #res = res[1]
    # println(string("res: ", res))

    for r = res
        bam = "/mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/bam.tsv"

        ref = "/mnt/raid64/Covid19_Gravida/cellranger/Homo_sapiens/genes/genes.sorted.gtf.gz"
    
        o = "/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_R/test1.pdf"
        o2 = "/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_R/test2.pdf"
        o3 = "/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_R/test3.pdf"
    
        lines = [1472189 - round(Int, i) for i = r.alpha_arr]
        lines = join(map(string, lines), ",")

        println(lines)
    
        run(`sashimiplot junc --gtf $ref --bam $bam --sj 1000 --junc 1:1471149:1472708 --ie 1,1  --ps RF --ssm R1 --fileout $o --trackline $lines --focus 1471469-1472189`)
    end
end

# test()

