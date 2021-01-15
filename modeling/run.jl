
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
            default = 0
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
    using_R::Bool=false, verbose::Bool=false,
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
    res = @showprogress "Reading... " pmap(beds) do b
        temp = Extract.get_record_from_bam(
            bam, b, 
            distance=distance, 
            min_reads=min_reads
        )
        return temp
    end

    # extract reads and drop duplicates
    temp_data = Vector{Extract.ExtractedData}()
    for t = res
        for r = t
            push!(temp_data, r)
        end
    end

    error_log = string(output, "_error.log")
    if exists(Path(error_log))
        rm(error_log)
    end

    res = @showprogress "Computing... "  pmap(collect(values(temp_data))) do data
        return Extract.run(
            data; n_max_ats = n_max_ats, n_min_ats = n_min_ats,
            mu_f = mu_f, sigma_f = sigma_f, min_ws = min_ws, 
            max_beta = max_beta, fixed_inference_flag = fixed_inference,
            single_end_mode = single_end, verbose = verbose, 
            using_R = using_R, error_log = error_log,
        )
    end

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
    max_beta=args["max-beta"], using_R=args["using-R"], 
    verbose=args["verbose"], process=args["process"], 
    n_max_ats=args["n-max-ats"], n_min_ats=args["n-min-ats"],
    fixed_inference=args["fixed-inference"], 
    single_end=args["single-end"], min_reads=args["min-reads"]
)