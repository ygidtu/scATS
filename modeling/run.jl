#!/usr/bin/env julia

using ArgParse
using BSON
using Distributed
using Memento
using FilePathsBase
using ProgressMeter


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Path to utr bed"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Prefix of output file"
            arg_type = String
            required = true
        "--using-R"
            help="How many process to use for R kernel, slow but more ATS sites. -1 means do not use R"
            arg_type = Int
            default = -1
        "--threads", "-t"
            help = "How many threads to use"
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
        "--chunk"
            help = "chunk size for multi threading, reduce memory usage."
            arg_type = Int64
            default = 1000
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
        "bam"
            help = "Path to bam file"
            required = true
            action = :store_arg
            nargs = '+'
    end

    return parse_args(s)
end

args = parse_commandline()

logger = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")
if args["verbose"]
    logger = Memento.config!("debug"; fmt="[{date} - {level} | {name}]: {msg} | {stacktrace}")
end

if args["using-R"] > 1
    addprocs(args["using-R"])
end

@everywhere begin
    using BGZFStreams
    using BioAlignments
    using FilePathsBase
    using ProgressMeter
    using XAM
    
    include(joinpath(@__DIR__, "src", "ATS.jl"))
    include(joinpath(@__DIR__, "src", "extract.jl"))
    include(joinpath(@__DIR__, "src", "genomic.jl"))

    function fit(
        b, bam::Vector;
        mu_f::Int64=300, min_ws::Float64=0.01,
        max_beta::Float64=50.0, sigma_f::Int64=50,
        n_max_ats::Int=5, n_min_ats::Int=1,
        fixed_inference::Bool=false, single_end::Bool=false,
        min_reads::Int=0, cage_mode::Bool = false,
        using_R::Bool=true, seed::Int=42,
        error_log::String = ""
    )
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
                using_R = using_R, error_log = error_log, seed=seed
            )
            return r
        end
        return ""
    end    
end



function normal_pipeline(
    input_file::String, bam::Vector, output::Union{String, PosixPath};
    mu_f::Int64=300, min_ws::Float64=0.01,
    max_beta::Float64=50.0, sigma_f::Int64=50,
    n_max_ats::Int=5, n_min_ats::Int=1,
    fixed_inference::Bool=false, single_end::Bool=false,
    min_reads::Int=0, cage_mode::Bool = false,
    using_R::Int=-1, seed::Int=42,
    chunk::Int = 1000
)   
    info(logger, "Loading UTR")
    beds = []  # Dict{String, Vector}()
    open(input_file, "r") do r
        seekend(r)
        fileSize = position(r)
        seekstart(r)
        p = Progress(fileSize)   # minimum update interval: .5 second
        # next!(p)
        while !eof(r)
            temp_bed = Genomic.new_bed(readline(r)) #, chromosomes=chromosomes)
            push!(beds, temp_bed)
            update!(p, position(r))

            if length(beds) >= 500
                break
            end
        end
        close(r)
    end

    info(logger, string("The number of UTRs: ", length(beds)))

    info(logger, "check index of bam")
    for b = bam
        bai = string(b, ".bai")
        if !exists(Path(bai))
            run(`samtools index $b`)
        end
    end
    
    info(logger, "start running")
    error_log = string(output, "_error.log")
    if exists(Path(error_log))
        rm(error_log)
    end

    header = ATSMIX.EMHeader()
    if using_R > 0
        res = @showprogress "Computing... " pmap(beds) do b
            return fit(
                b, bam, 
                n_max_ats=n_max_ats, n_min_ats=n_min_ats,
                mu_f = mu_f, sigma_f = sigma_f, min_ws = min_ws, 
                max_beta = max_beta, fixed_inference = fixed_inference,
                single_end = single_end, cage_mode=cage_mode,
                using_R = using_R > 0, error_log = error_log, seed=seed,
                min_reads=min_reads
            )
        end
        w = open(output, "w+")
        write(w, string(join(header, "\t"), "\n"))

        @showprogress 1 "Writing..." for r = res
            if r != ""
                write(w, string(r, "\n"))
            end
        end
        close(w)
    else
        open(output, "w+") do w
            
            write(w, string(join(header, "\t"), "\n"))

            p = Progress(length(beds), 1, "Computing...")

            for i = 1:chunk:length(beds)
                Threads.@threads for j = i:min(i + chunk - 1, length(beds))
                    r = fit(
                        beds[j], bam, 
                        n_max_ats=n_max_ats, n_min_ats=n_min_ats,
                        mu_f = mu_f, sigma_f = sigma_f, min_ws = min_ws, 
                        max_beta = max_beta, fixed_inference = fixed_inference,
                        single_end = single_end, cage_mode=cage_mode,
                        using_R = using_R > 0, error_log = error_log, seed=seed,
                        min_reads=min_reads
                    )

                    if r != ""
                        write(w, string(join([r[x] for x = header], "\t"), "\n"))
                    end

                    r = nothing

                    GC.safepoint()
                    next!(p)
                end
                if i > chunk
                    Base.GC.gc()
                end
            end
            close(w)
        end
    end
end


function identify_pipeline(
    input_file::String, bam::Vector, output::Union{String, PosixPath};
    mu_f::Int64=300, min_ws::Float64=0.01,
    max_beta::Float64=50.0, sigma_f::Int64=50,
    n_max_ats::Int=5, n_min_ats::Int=1,
    fixed_inference::Bool=false, single_end::Bool=false,
    min_reads::Int=0, cage_mode::Bool = false,
    expand::String="200,800", 
    using_R::Bool=true, seed::Int=42
)
    gtf = Genomic.load_GTF(input_file)

    res = @showprogress "Computing... " pmap(collect(keys(gtf["gene"]))) do g_name
        for t_name = gtf["gene"][g_name]
            transcript = gtf["transcript"][t_name]
            data = Extract.get_record_from_bam_transcript(
                bam, 
                transcript, 
                get(gtf["exon"], t_name, []),
                expand = expand
            )

            if isnothing(data) || length(data.st_arr) == 0
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
    
                return r
            end
            return ""
        end
    end
    w = open(output, "w+")
    write(w, "gene_id\ttranscript_id\tutr\tst_arr\ten_arr\tws\talpha_arr\tbeta_arr\tlb_arr\tlabel\tbic\n")
    for r = res
        write(w, string(g_name, "\t", t_name, "\t", r, "\n"))
    end
    close(w)
end


function main(
    input_file::String, bam::Vector, output::String;
    mu_f::Int64=300, min_ws::Float64=0.01,
    max_beta::Float64=50.0, sigma_f::Int64=50,
    using_R::Int=-1, chunk::Int=1000,
    n_max_ats::Int=5, n_min_ats::Int=1,
    fixed_inference::Bool=false, single_end::Bool=false,
    min_reads::Int=0, cage_mode::Bool = false,
    identify::Bool = false, expand::String="200,1000",
    debug::Bool = false, seed::Int=42
)
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

    if using_R == 0
        using_R = 1
    end

    if identify
        if debug
            return test1()
        end
        identify_pipeline(
            input_file, bam, output, 
            mu_f=mu_f, min_ws=min_ws,
            max_beta=max_beta, sigma_f=sigma_f,
            using_R=using_R > 0,
            n_max_ats=n_max_ats, n_min_ats=n_min_ats,
            fixed_inference=fixed_inference, single_end=single_end,
            min_reads=min_reads, cage_mode=cage_mode,
            expand=expand, seed=seed
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
            using_R=using_R,
            n_max_ats=n_max_ats, n_min_ats=n_min_ats,
            fixed_inference=fixed_inference, single_end=single_end,
            min_reads=min_reads, cage_mode=cage_mode,
            seed=seed, chunk=chunk
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

    function generate_intron(sites, exon_coord, exon_range)
        relative_range = [exon_coord[x] for x = exon_range]

        split_sites = []
   
        for j = 1:2:length(exon_range)
            if relative_range[j] <= sites[1] < sites[2] <= relative_range[j + 1]
                push!(split_sites, exon_range[j] + sites[1] - relative_range[j])
                push!(split_sites, exon_range[j] + sites[2] - relative_range[j])
            elseif relative_range[j] <= sites[1] <= relative_range[j + 1] < sites[2]
                push!(split_sites, exon_range[j] + sites[1] - relative_range[j])
                push!(split_sites, exon_range[j + 1])
            elseif sites[1] < relative_range[j] <= sites[2] <= relative_range[j + 1]
                push!(split_sites, exon_range[j])
                push!(split_sites, exon_range[j] + sites[2] - relative_range[j])
            end


            if relative_range[j] <= sites[3] < sites[4] <= relative_range[j + 1]
                push!(split_sites, exon_range[j] + sites[3] - relative_range[j])
                push!(split_sites, exon_range[j] + sites[4] - relative_range[j])
            elseif relative_range[j] <= sites[3] <= relative_range[j + 1] < sites[4]
                push!(split_sites, exon_range[j] + sites[1] - relative_range[j])
                push!(split_sites, exon_range[j + 1])
            elseif sites[3] < relative_range[j] <= sites[4] <= relative_range[j + 1]
                push!(split_sites, exon_range[j])
                push!(split_sites, exon_range[j] + sites[4] - relative_range[j])
            end
        end
 
        return split_sites
    end

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
                        break
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
    
    for seed = 1:1
        gtf = Genomic.load_GTF(joinpath(@__DIR__, "..", string("tests/", seed, ".gtf")))

        real_sites = Vector()
        junc_sites = Vector()

        #=
        t_name = "T1"
        exons = gtf["exon"][t_name]
        # convert genomic site of exons to relative pos in transcript
        exon_coord = Dict{Int, Int}()
        exon_range = Vector()
        for e = exons
            for i = e.Start:e.End
                exon_coord[i] = exons[1].Start + length(exon_coord) + 1 - gtf["transcript"][t_name].Start
            end
            append!(exon_range, [e.Start, e.End])
        end
        println(exon_range)
        println([exon_coord[x] for x = exon_range])

        is_read1_junc = true # || line[6] == "TRUE"
        is_read2_junc = false

        line = [3002, 3141, 4391, 4537]
        # println(line)
        r1_site = [line[1], line[2]]
        r2_site = [line[3], line[4]]

        site_relative = [
            convert_absolute_relative(exon_coord, exon_range, r1_site, is_read1_junc, is_read2 = false)...,
            convert_absolute_relative(exon_coord, exon_range, r2_site, is_read2_junc, is_read2 = true)...,
        ]

        println(site_relative)
        println(generate_intron(site_relative, exon_coord, exon_range))
        =#

        open(joinpath(@__DIR__, "..", string("tests/sites_", seed, ".txt"))) do r
            while !eof(r)
                line = readline(r)
                line = split(strip(line), ",")
                
                t_name = line[7]

                exons = gtf["exon"][t_name]
                # convert genomic site of exons to relative pos in transcript
                exon_coord = Dict{Int, Int}()
                exon_range = Vector()
                for e = exons
                    for i = e.Start:e.End
                        exon_coord[i] = exons[1].Start + length(exon_coord) + 1 - gtf["transcript"][t_name].Start
                    end
                    append!(exon_range, [e.Start, e.End])
                end

                is_read1_junc = line[5] == "TRUE" # || line[6] == "TRUE"
                is_read2_junc = line[6] == "TRUE"

                line = parse.(Int, line[1:4])
                # println(line)
                r1_site = [line[1], line[2]]
                r2_site = [line[3], line[4]]

                site_relative = [
                    convert_absolute_relative(exon_coord, exon_range, r1_site, is_read1_junc, is_read2 = false)...,
                    convert_absolute_relative(exon_coord, exon_range, r2_site, is_read2_junc, is_read2 = true)...,
                ]
                append!(line, [is_read1_junc ? 1 : 0, is_read2_junc ? 1 : 0])
                push!(real_sites, line)
                push!(junc_sites, generate_intron(site_relative, exon_coord, exon_range))
            end
            close(r)
        end

        for (t_name, transcript) = gtf["transcript"]
            exons = gtf["exon"][t_name]
            # convert genomic site of exons to relative pos in transcript
            exon_coord = Dict{Int, Int}()
            exon_range = Vector()
            for e = exons
                for i = e.Start:e.End
                    exon_coord[i] = exons[1].Start + length(exon_coord) + 1 - gtf["transcript"][t_name].Start
                end
                append!(exon_range, [e.Start, e.End])
            end


            st_arr = Vector()
            en_arr = Vector()
            r1_st = Vector()
            r1_en = Vector()
            r2_st = Vector()
            r2_en = Vector()
            relative = Vector()

            # r1_site, r2_site = [2336, 2442], [2701, 2792]

            # is_read1_junc = false
            # is_read2_junc = true
            # site_relative = [
            #     convert_absolute_relative(exon_coord, exon_range, r1_site, is_read1_junc, is_read2 = false)...,
            #     convert_absolute_relative(exon_coord, exon_range, r2_site, is_read2_junc, is_read2 = true)...,
            # ]
            # println(site_relative)

            for (sites, junc) = zip(real_sites, junc_sites)
                debug_site = [
                    # [2336, 2792],
                    # [1934, 3947],
                    # [2554, 3957],
                    # [1952, 2723],
                    # [2970, 4467],
                    # [3002, 4537],
                    # [2268, 4190],
                    [2614, 4268]
                ]

                pass = false
                failed_on_border = true
                pass1 = false
                pass2 = false

                # println(junc_sites)
                relative_exons = [exon_coord[x] for x = exon_range]
                if site_relative[2] in relative_exons || site_relative[3] in relative_exons
                    failed_on_border = abs(site_relative[3] - site_relative[2]) <= 1
                else
                    failed_on_border = false
                end
                
                
                for j = 1:2:length(relative_exons)
                    if is_read1_junc
                        if relative_exons[j] <= site_relative[1] <= relative_exons[j + 1] < site_relative[2]
                            pass1 = true
                        end
                    elseif relative_exons[j] <= site_relative[1] < site_relative[2] <= relative_exons[j + 1]
                        pass1 = true
                    end

                    if is_read2_junc
                        if site_relative[3] < relative_exons[j] <= site_relative[4] <= relative_exons[j + 1]
                            pass2 = true
                        end
                    elseif relative_exons[j] <= site_relative[3] < site_relative[4] <= relative_exons[j + 1]
                        pass2 = true
                    end
                    
                    if pass1 && pass2
                        break
                    end
                end

                
                pass = pass1 && pass2 && !failed_on_border

                # println(sites)
                relative = [
                    convert_absolute_relative(exon_coord, exon_range, sites[1:2], sites[5] == 1, is_read2 = false)...,
                    convert_absolute_relative(exon_coord, exon_range, sites[3:4], sites[6] == 1, is_read2 = true)...
                ]

                if any([sites[1] == x[1] && sites[4] == x[2] for x = debug_site])
                    println(t_name)
                    println(sites)
                    println(relative)
                    println(junc)
                    println([haskey(exon_coord, x) for x = junc])
                end

                if !any(isnan.(relative)) && all([haskey(exon_coord, x) for x = junc])
                    push!(st_arr, get(exon_coord, sites[1], 0))
                    push!(en_arr, get(exon_coord, sites[4], 0))

                    push!(r1_st, sites[1])
                    push!(r1_en, sites[2])
                    push!(r2_st, sites[3])
                    push!(r2_en, sites[4])

                    # push!(relative, join(map(string, site_relative), "|"))
                end
            end

            df = DataFrame(st_arr = st_arr, en_arr = en_arr, r1_st = r1_st, r1_en = r1_en, r2_st = r2_st, r2_en = r2_en)

            CSV.write(joinpath(@__DIR__, "..", "tests", string(t_name, "_", seed, ".csv")), df)
        end
    end
end



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
    identify = args["identify"],
    seed=args["seed"], chunk=args["chunk"],
    debug=get(args, "debug", false)
)
