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
    
    using Gadfly
    using Cairo
    using DataFrames
    using CSV

    include(joinpath(@__DIR__, "src", "ATS.jl"))
    include(joinpath(@__DIR__, "src", "extract.jl"))
    include(joinpath(@__DIR__, "src", "genomic.jl"))
end


function normal_pipeline(
    input_file::String, bam::String, output::Union{String, PosixPath};
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


function identify_pipeline(
    input_file::String, bam::String, output::Union{String, PosixPath};
    mu_f::Int64=300, min_ws::Float64=0.01,
    max_beta::Float64=50.0, sigma_f::Int64=50,
    using_R::Bool=true, verbose::Bool=false,
    n_max_ats::Int=5, n_min_ats::Int=1,
    fixed_inference::Bool=false, single_end::Bool=false,
    min_reads::Int=0, cage_mode::Bool = false,
    expand::String="200,800"
)
    gtf = Genomic.load_GTF(input_file)

    w = open(output, "w+")
    @showprogress "Computing... " pmap(collect(keys(gtf["gene"]))) do g_name
        for t_name = gtf["gene"][g_name]
            transcript = gtf["transcript"][t_name]
            data = Extract.get_record_from_bam_transcript(
                bam, 
                transcript, 
                get(gtf["exon"], t_name, []),
                expand = expand
            )

            if length(data.st_arr) == 0
                continue
            end

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
                    debug = false
                )
    
                if r != ""
                    write(w, "gene_id\ttranscript_id\tutr\tst_arr\ten_arr\tws\talpha_arr\tbeta_arr\tlb_arr\tlabel\tbic\n")
                    write(w, string(g_name, "\t", t_name, "\t", r, "\n"))
                end
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
    identify::Bool = false, expand::String="200,1000",
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
        identify_pipeline(
            input_file, bam, output, 
            mu_f=mu_f, min_ws=min_ws,
            max_beta=max_beta, sigma_f=sigma_f,
            using_R=using_R, verbose=verbose,
            n_max_ats=n_max_ats, n_min_ats=n_min_ats,
            fixed_inference=fixed_inference, single_end=single_end,
            min_reads=min_reads, cage_mode=cage_mode,
            expand=expand
        )
    else
        if debug
            # return test()
            return test_draw_density()
        end
        normal_pipeline(
            input_file, bam, output, 
            mu_f=mu_f, min_ws=min_ws,
            max_beta=max_beta, sigma_f=sigma_f,
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
            bam = "/mnt/raid61/Personal_data/zhangyiming/code/afe/tests/bam.tsv"

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


function test1(mu_f::Int64=300, min_ws::Float64=0.01,
    max_beta::Float64=50.0, sigma_f::Int64=70,
    using_R::Bool=false, cage_mode=false,
    n_max_ats::Int=5, n_min_ats::Int=2,
    fixed_inference::Bool=false, single_end::Bool=false,
    min_reads::Int=0)

    function convert_absolute_relative(exon_coord::Dict, exon_range::Vector, sites::Vector, is_reads_junc::Bool; is_read2::Bool = false)::Vector
        site_relative = []
        if is_reads_junc
            if is_read2
                s = NaN
                for i = 1:2:length(exon_range)
                    if exon_range[i] <= sites[2] <= exon_range[i + 1]
                        s = exon_coord[exon_range[i]] - (exon_range[i] - sites[1])
                        break
                    end
                end
                push!(site_relative, s)
                push!(site_relative, get(exon_coord, sites[2], NaN))

            else
                push!(site_relative, get(exon_coord, sites[1], NaN))

                s = NaN
                for i = 1:2:length(exon_range)
                    if exon_range[i] <= sites[1] <= exon_range[i + 1]
                        s = exon_coord[exon_range[i + 1]] + sites[2] - exon_range[i + 1]
                    end
                end
                push!(site_relative, s)
            end
        else
            for j = sites
                push!(site_relative, get(exon_coord, j, NaN))
            end
        end
 
        return site_relative
    end

    gtf = Genomic.load_GTF(joinpath(@__DIR__, "..", "tests/sim.gtf"))

    t_name = "T1"
    transcript = gtf["transcript"][t_name]

    for (t_name, transcript) = gtf["transcript"]

        utr = transcript.Start
        exons = gtf["exon"][t_name]
        # convert genomic site of exons to relative pos in transcript
        exon_coord = Dict{Int, Int}()
        exon_range = Vector()
        for e = exons
            for i = e.Start:e.End
                exon_coord[i] = exons[1].Start + length(exon_coord) + 1 - utr
            end
            append!(exon_range, [e.Start, e.End])
        end


        st_arr = Vector()
        en_arr = Vector()
        real_st = Vector()
        real_en = Vector()

        # r1_site, r2_site = [2336, 2442], [2701, 2792]

        # is_read1_junc = false
        # is_read2_junc = true
        # site_relative = [
        #     convert_absolute_relative(exon_coord, exon_range, r1_site, is_read1_junc, is_read2 = false)...,
        #     convert_absolute_relative(exon_coord, exon_range, r2_site, is_read2_junc, is_read2 = true)...,
        # ]
        # println(site_relative)
        open(joinpath(@__DIR__, "..", "tests/sites.txt")) do r
            while !eof(r)
                line = readline(r)
                line = split(strip(line), ",")
                is_read1_junc = line[5] == "TRUE" # || line[6] == "TRUE"
                is_read2_junc = line[6] == "TRUE"
                line = parse.(Int, line[1:4])
                # println(line)
                r1_site = [line[1], line[2]]
                r2_site = [line[3], line[4]]

                #=
                line = [2261,2357,4201,4295]
                r1_site = [2261, 2357]
                r2_site = [4201,4295]
                =#

                start_site, end_site = r1_site[1], r2_site[2]
                
                site_relative = [
                    convert_absolute_relative(exon_coord, exon_range, r1_site, is_read1_junc, is_read2 = false)...,
                    convert_absolute_relative(exon_coord, exon_range, r2_site, is_read2_junc, is_read2 = true)...,
                ]

                debug_site = [
                    # [2336, 2792],
                    # [1934, 3947],
                    # [2554, 3957],
                    [2257, 4045],
                    [2255, 4002]
                ]

                if !any(isnan.(site_relative))
                    pass = false

                    if is_read1_junc || is_read2_junc # && (site_relative[3] - site_relative[2] > 1 && any([abs(exon_coord[x] - site_relative[2]) <= 1 for x = exon_range]))
                        for i = 1:2:length(exon_range)
                            if exon_coord[exon_range[i]] <= site_relative[2] < exon_coord[exon_range[i + 1]]
                                if exon_coord[exon_range[i]] < site_relative[3] <= exon_coord[exon_range[i+1]]
                                    
                                    # if reads is just locate on the edge of exons and junctions, then R1 and R2 needs to have at least 2 bp distance
                                    relative_exons = [exon_coord[x] for x = exon_range]
                                    if site_relative[2] in relative_exons || site_relative[3] in relative_exons
                                        pass = site_relative[3] - site_relative[2] > 1
                                    else
                                        pass = true
                                    end
                                    break
                                end
                            end
                        end
                    end
                    
                    # if r1_site[1] == 2257 && r2_site[2] == 4045 || r1_site[1] == 2083 && r2_site[2] == 4050 || r1_site[1] == 2239 && r2_site[2] == 4361
                    if any([r1_site[1] == x[1] && r2_site[2] == x[2] for x = debug_site])
                        println(string(t_name, ": ", is_read1_junc, "-", is_read2_junc, "; pass = ", pass, "; sites: ", site_relative))
                        println(string("r1: ", r1_site, "; r2: ", r2_site))
                        println(string("!(is_read1_junc || is_read2_junc) = ", !(is_read1_junc || is_read2_junc), " all(site_relative .> 0) == ", all(site_relative .> 0)))
                        println(string("exon_range: ", [exon_coord[x] for x = exon_range]))
                    end

                    if !(is_read1_junc || is_read2_junc) || pass
                        push!(st_arr, site_relative[1])
                        push!(real_st, start_site)

                        push!(en_arr, site_relative[4])
                        push!(real_en, end_site)
                        continue
                    end
                else
                    if any([r1_site[1] == x[1] && r2_site[2] == x[2] for x = debug_site])
                        println(string(t_name, ": ", is_read1_junc, "-", is_read2_junc, "; sites: ", site_relative))
                        println(string("r1: ", r1_site, "; r2: ", r2_site))
                        println(string("!(is_read1_junc || is_read2_junc) = ", !(is_read1_junc || is_read2_junc), " all(site_relative .> 0) == ", all(site_relative .> 0)))
                        println(string("exon_range: ", [exon_coord[x] for x = exon_range]))
                    end
                end
   
                push!(st_arr, 0)
                push!(real_st, start_site)

                push!(en_arr, 0)
                push!(real_en, end_site)
            end
            close(r)
        end

        df = DataFrame(st_arr = st_arr, en_arr = en_arr, real_st = real_st, real_en = real_en)

        CSV.write(joinpath(@__DIR__, "..", "tests", string(t_name, ".csv")), df)
    end
end


    # bam = "/mnt/raid64/ATS/Personal/zhangyiming/bams/NHC2.bam"
    # ref = "/mnt/raid64/Covid19_Gravida/cellranger/Homo_sapiens/genes/genes.sorted.gtf.gz"  
    # o = "/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_jl/"
    # gtf = Genomic.load_GTF(joinpath(@__DIR__, "..", "tests/genes.gtf"))
    # seed = 42
    # res = [[], []]
    # for (t_name, transcript) = gtf["transcript"]
    #     println(transcript)
    #     println([e.Start for e = get(gtf["exon"], t_name, [])])
    #     data = Extract.get_record_from_bam_transcript(bam, transcript, get(gtf["exon"], t_name, []))

    #     if !isnothing(data)
    #         # r = ATSMIX.fit(
    #         #     data.utr,
    #         #     n_max_ats, n_min_ats,
    #         #     data.st_arr, data.en_arr,
    #         #     mu_f = mu_f, sigma_f = sigma_f, min_ws = min_ws, 
    #         #     max_beta = max_beta, fixed_inference_flag = fixed_inference,
    #         #     single_end_mode = single_end,
    #         #     using_R = using_R, error_log = nothing, seed=seed,
    #         #     cage_mode=cage_mode, exon_coord = data.exon_coord,
    #         #     debug = true
    #         # )

    #         # if r != ""
    #         #     junc = string(data.utr.Chrom, ":", data.utr.Start - 1000, ":", data.utr.End + 1000)
    #         #     lines = join(map(string, r.absolute_arr), ",")
    #         #     focus = string(data.utr.Start, "-", data.utr.End)

    #         #     tsv = "/mnt/raid61/Personal_data/zhangyiming/code/afe/tests/bam.tsv"

    #         #     run(`sashimiplot junc --gtf $ref --bam $tsv --sj 1000 --junc $junc --ie 1,1  --ps RF --ssm R1 --fileout $o/$t_name.pdf --trackline $lines --focus $focus`)  # 1212795-1214738
    #         # end

    #         append!(res[2], abs.(data.st_arr .- data.en_arr))
    #         append!(res[1], [Symbol(t_name) for _ = 1:length(data.st_arr)])
    #     end
    # end

    # df = DataFrame(Transcript = res[1], Len = res[2])

    # @df df density(:Len, group = (:Transcript), legend = :topright)
    # savefig(string(o, "/ENST_density.pdf"))

function test_draw_density()
    bam = "/mnt/raid64/ATS/Personal/zhangyiming/bams/NHC2.bam"
    ref = "/mnt/raid64/Covid19_Gravida/cellranger/Homo_sapiens/genes/genes.sorted.gtf.gz"  
    o = "/mnt/raid64/ATS/Personal/zhangyiming/CG/NHC2_jl_all/test_jl/genes_density"
    # gtf = Genomic.load_GTF(joinpath(@__DIR__, "..", "tests/genes.gtf"))
    gtf = Genomic.load_GTF("/mnt/raid64/Covid19_Gravida/cellranger/Homo_sapiens/genes/genes.gtf")
    seed = 42
    
    # println(collect(keys(gtf["gene"])))
    @showprogress pmap(collect(keys(gtf["gene"]))) do g_name
    # for g_name = collect(keys(gtf["gene"]))
        t_names = gtf["gene"][g_name]

        res = [[], []]
        for t_name = t_names
            transcript = gtf["transcript"][t_name]
            data = Extract.get_record_from_bam_transcript(bam, transcript, get(gtf["exon"], t_name, []), expand="200,500")

            if !isnothing(data)
                append!(res[2], abs.(data.st_arr .- data.en_arr))
                append!(res[1], [Symbol(t_name) for _ = 1:length(data.st_arr)])
            end
        end

        df = DataFrame(Transcript = res[1], Len = res[2])

        if size(df)[1] > 0
            # @df df density(:Len, group = (:Transcript), legend = :topright)
            # savefig(joinpath(o, string(g_name, ".pdf")))
            p = plot(df, x=:Len, color=:Transcript,
                Geom.density, 
                Guide.ylabel("Density"), 
                Theme(alphas=[0.6]),
                Guide.title(g_name),
                Scale.x_continuous(minvalue=-0, maxvalue=min(1500, maximum(df[:Len])))
            )
            draw(PDF(joinpath(o, string(g_name, ".pdf")), 5inch, 3inch), p)
        end
    end
end



main(
    args["input"], args["bam"], args["output"], 
    mu_f = args["mu-f"], min_ws = args["min-ws"], 
    sigma_f = args["sigma-f"], max_beta=args["max-beta"], 
    using_R=args["using-R"], cage_mode=args["cage-mode"],
    n_max_ats=args["n-max-ats"], n_min_ats=args["n-min-ats"],
    fixed_inference=args["fixed-inference"], single_end=args["single-end"], 
    min_reads=args["min-reads"], expand=args["expand"],
    verbose=args["verbose"], identify = args["identify"],
    debug=get(args, "debug", false)
)

