#!/usr/bin/env julia

using ArgParse
using BioGenerics
using BSON
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
            default = 0.1
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
            default = 5
        "--expand"
            help = "the upstream and downstream UTR region around start site of each transcript"
            arg_type = String
            default = "200,1000"
        "--single-end"
            help="whether this is sinle-end sequencing"
            action = :store_true
        "--fixed-inference"
            help="inference with fixed parameters"
            action = :store_true
        "--cage-mode"
            help="whether to run in cage mode to identify ATS, only work with PE data"
            action = :store_true
        "--identify"
            help="whether to identify the host transcript of this ATS"
            action = :store_true
        "--verbose"
            help="whether to dispaly more detailed information"
            action = :store_true
        "--debug"
            help="whether to run in debug mode"
            action = :store_true
    end

    return parse_args(s)
end

args = parse_commandline()

if args["process"] > 1
    addprocs(args["process"])
end

@everywhere begin
    using FilePathsBase
    using Memento
    using ProgressMeter
    
    include(joinpath(@__DIR__, "src", "ATS.jl"))
    include(joinpath(@__DIR__, "src", "extract.jl"))
    include(joinpath(@__DIR__, "src", "genomic.jl"))
end


function normal_pipeline(
    input_file::String, bam::String, output::String;
    mu_f::Int64=300, min_ws::Float64=0.01,
    max_beta::Float64=50.0, sigma_f::Int64=50,
    using_R::Bool=true, verbose::Bool=false,
    n_max_ats::Int=5, n_min_ats::Int=1,
    fixed_inference::Bool=false, single_end::Bool=false,
    min_reads::Int=0, cage_mode::Bool = false,
)
    beds = []  # Dict{String, Vector}()
    r = open(input_file, "r")
    seekend(r)
    fileSize = position(r)
    seekstart(r)
    p = Progress(fileSize, 1)   # minimum update interval: 1 second
    while !eof(r)
        temp_bed = Genomic.new_bed(readline(r)) #, chromosomes=chromosomes)

        push!(beds, temp_bed)

        # if length(beds) > 200
        #     break
        # end

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
    
    w = open(output, "w+")
    write(w, "utr\tst_arr\ten_arr\tws\talpha_arr\tbeta_arr\tlb_arr\tlabel\tbic\n")
    @showprogress "Computing... " pmap(beds) do b
        data = Extract.get_record_from_bam(
            bam, b, 
            min_reads=min_reads, 
            cage_mode=cage_mode, 
            single_end_mode=single_end
        )

        if !isnothing(data)
            r = ATSMIX.fit(
                data.utr, n_max_ats, n_min_ats, data.st_arr, data.en_arr,
                mu_f = mu_f, sigma_f = sigma_f, min_ws = min_ws, 
                max_beta = max_beta, fixed_inference_flag = fixed_inference,
                single_end_mode = single_end, cage_mode=cage_mode,
                using_R = using_R, error_log = error_log
            )
            if r != ""
                write(w, string(r, "\n"))
            end
        end
    end
    close(w)
end



function main(
    input_file::String, bam::String, output::String;
    mu_f::Int64=300, min_ws::Float64=0.01,
    max_beta::Float64=50.0, sigma_f::Int64=50,
    using_R::Bool=true, verbose::Bool=false,
    n_max_ats::Int=5, n_min_ats::Int=1,
    fixed_inference::Bool=false, single_end::Bool=false,
    min_reads::Int=0, cage_mode::Bool = false,
    identify::Bool = false,
    debug::Bool = false
)

    logger = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")
    if verbose
        logger = Memento.config!("debug"; fmt="[{date} - {level} | {name}]: {msg} | {stacktrace}")
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

    if verbose
        Extract.setLevel("debug")
    end

    if identify
        if debug
            return test1()
        end
    else
        if debug
            return test()
        end
        normal_pipeline(
            input_file, bam, output, 
            mu_f=mu_f, min_ws=min_ws,
            max_beta=meta_beta, sigma_f=sigma_f,
            using_R=using_R, verbose=verbose,
            n_max_ats=n_max_ats, n_min_ats=n_min_ats,
            fixed_inference=fixed_inference, single_end=single_end,
            min_reads=min_reads, cage_mode=cage_mode,
        )
    end
end



function test(mu_f::Int64=300, min_ws::Float64=0.01,
    max_beta::Float64=50.0, sigma_f::Int64=70,
    using_R::Bool=true, cage_mode=false,
    n_max_ats::Int=5, n_min_ats::Int=2,
    fixed_inference::Bool=false, single_end::Bool=false,
    min_reads::Int=0)
    println("run test func")
    bam = "/mnt/raid64/ATS/Personal/zhangyiming/bams/NHC2.bam"

    bed = Genomic.new_bed("1\t1212595\t1214738\t.\t.\t-")
    # "1\t6205375\t6206201\t.\t.\t+"
    # 1\t2556266\t2556833\t.\t.\t+
    # 1\t1471469\t1472189\t.\t.\t+
    
    focus = Vector()
    open("/mnt/raid64/ATS/Personal/zhangyiming/cage/CG_R1.bed") do r
        while !eof(r)
            line = readline(r)
            b = Genomic.new_bed(line)

            if b.Chrom == "1" && b.Start < 1214738 && b.End > 1212595
                push!(focus, string(b.Start, "-", b.End))
            end
        end
        close(r)
    end

    temp = Extract.get_record_from_bam( bam, bed, min_reads=min_reads, cage_mode=cage_mode)

    if isnothing(temp)
        return
    end

    # pmap(1:20) do seed
    for seed = 1:4
        r = ATSMIX.fit(
            temp.utr,
            n_max_ats, n_min_ats,
            temp.st_arr, temp.en_arr,
            mu_f = mu_f, sigma_f = sigma_f, min_ws = min_ws, 
            max_beta = max_beta, fixed_inference_flag = fixed_inference,
            single_end_mode = single_end,
            using_R = using_R, error_log = nothing, seed=seed,
            cage_mode=cage_mode, debug = true
        )

        if !isnothing(r.alpha_arr) && length(r.alpha_arr)  > 0
            println(string("seed: ", seed))
            bam = "/mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/bam.tsv"

            ref = "/mnt/raid64/Covid19_Gravida/cellranger/Homo_sapiens/genes/genes.sorted.gtf.gz"
        
            o = string("/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_jl/test1_", seed, ".pdf")
            o2 = "/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_jl/test2.pdf"
            o3 = "/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_jl/test3.pdf"
        
            lines = [1212795, [round(Int, i) for i = r.absolute_arr]...]
            lines = join(map(string, lines), ",")

            # dots = join(map(string, temp[1].real_st), ",")
            # dots2 =  join(map(string, temp[1].real_en), ",")

            # focus = join(focus, ",")
            # print(dots)
            run(`sashimiplot junc --gtf $ref --bam $bam --sj 10 --junc 1:1210795:1214738 --ie 1,1  --ps RF --ssm R1 --fileout $o2 --focus 1212795-1214738 --trackline $lines`)  #   --dots $dots
            break
        end

        
    end

end


using Plots
using StatsPlots
using DataFrames

function test1(mu_f::Int64=300, min_ws::Float64=0.01,
    max_beta::Float64=50.0, sigma_f::Int64=70,
    using_R::Bool=false, cage_mode=false,
    n_max_ats::Int=5, n_min_ats::Int=2,
    fixed_inference::Bool=false, single_end::Bool=false,
    min_reads::Int=0)
    bam = "/mnt/raid64/ATS/Personal/zhangyiming/bams/NHC2.bam"

    gtf = Genomic.load_GTF(joinpath(@__DIR__, "genes.gtf"))
    ref = "/mnt/raid64/Covid19_Gravida/cellranger/Homo_sapiens/genes/genes.sorted.gtf.gz"
        
    o = "/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_jl/"

    seed = 42
    res = [[], []]
    for (t_name, transcript) = gtf["transcript"]
        println(transcript)
        println([e.Start for e = get(gtf["exon"], t_name, [])])
        data = Extract.get_record_from_bam_transcript(bam, transcript, get(gtf["exon"], t_name, []))

        if !isnothing(data)
            r = ATSMIX.fit(
                data.utr,
                n_max_ats, n_min_ats,
                data.st_arr, data.en_arr,
                mu_f = mu_f, sigma_f = sigma_f, min_ws = min_ws, 
                max_beta = max_beta, fixed_inference_flag = fixed_inference,
                single_end_mode = single_end,
                using_R = using_R, error_log = nothing, seed=seed,
                cage_mode=cage_mode, exon_coord = data.exon_coord,
                debug = true
            )

            # if r != ""
            #     junc = string(data.utr.Chrom, ":", data.utr.Start - 1000, ":", data.utr.End + 1000)
            #     lines = join(map(string, r.absolute_arr), ",")
            #     focus = string(data.utr.Start, "-", data.utr.End)

            #     tsv = "/mnt/raid61/Personal_data/zhangyiming/code/afe/modeling/bam.tsv"

            #     run(`sashimiplot junc --gtf $ref --bam $tsv --sj 1000 --junc $junc --ie 1,1  --ps RF --ssm R1 --fileout $o/$t_name.pdf --trackline $lines --focus $focus`)  # 1212795-1214738
            # end

            append!(res[2], abs.(data.st_arr .- data.en_arr))
            append!(res[1], [Symbol(t_name) for _ = 1:length(data.st_arr)])
        end
    end

    df = DataFrame(Transcript = res[1], Len = res[2])

    @df df density(:Len, group = (:Transcript), legend = :topright)
    savefig(string(o, "/ENST_density.pdf"))
end

main(
    args["input"], args["bam"], args["output"], 
    mu_f = args["mu-f"], min_ws = args["min-ws"], 
    sigma_f = args["sigma-f"], max_beta=args["max-beta"], 
    using_R=args["using-R"], cage_mode=args["cage-mode"],
    n_max_ats=args["n-max-ats"], n_min_ats=args["n-min-ats"],
    fixed_inference=args["fixed-inference"], single_end=args["single-end"], 
    min_reads=args["min-reads"],
    verbose=args["verbose"], identify = args["identify"],
    debug=get(args, "debug", false)
)

