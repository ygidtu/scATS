const jobs = RemoteChannel(()->Channel{Dict}(32));
const results = RemoteChannel(()->Channel{Union{Nothing, Dict}}(32));


@everywhere begin
    using ProgressMeter

    include(joinpath(@__DIR__, "src", "bam.jl"))
    include(joinpath(@__DIR__, "src", "core.jl"))
    include(joinpath(@__DIR__, "src", "genomic.jl"))


    function do_work(jobs, results)
        while true
            args = take!(jobs)
            st_arr = Bam.load_reads(args["bams"], args["utr"])

            if length(st_arr) < args["min-reads"]
                put!(results, nothing)
                continue
            end

            param = ATSModel.init(
                st_arr, 
                args["utr"],
                n_max_ats=args["n-max-ats"], 
                n_min_ats=args["n-min-ats"], 
                min_ws=args["min-ws"], 
                max_unif_ws=args["max-unif-ws"], 
                fixed_inference_flag=args["fixed-inference"],
                step_size=args["step-size"],
                nround=args["nround"],
                max_beta=args["max-beta"],
                seed=args["seed"]
            )

            try
                res = ATSModel.run(param)
                res = ATSModel.toDict(res)
                put!(results, res)
            catch err
                put!(results, nothing)
            end
        end     
    end
end


function make_jobs(utrs::Vector, args::Dict, bams::Vector)
    for i = utrs
        args["utr"] = i
        args["bams"] = bams
        put!(jobs, args)
    end
end


function infer(args::Dict; logger)
    bams = Bam.prepare_bam_list(args["bam"])

    info(logger, "Load UTR from GTF")
    utrs = Genomic.load_UTRs(args["gtf"], utr_length = args["utr-length"])

    output = absolute(Path(args["output"]))
    out_dir = parent(output)
    try
        if !exists(out_dir)
            mkdir(Path(out_dir), recursive=true)
        end
    catch e
        error(logger, Formatting.format(FormatExpr("Error while create {}: {}"), output, e))
        exit(1)
    end

    # Counting by multi-threads
    p = Progress(length(utrs), 1, "Computing...")
    w = open(output, "w+")

    # write header
    header = ATSModel.getHeader()
    write(w, string(join(header, "\t"), "\n"))
    flush(w)

    # multi-processing
    errormonitor(@async make_jobs(utrs, args, bams))

    for p = workers()
        remote_do(do_work, p, jobs, results)
    end

    n = length(utrs)
    @elapsed while n > 0
        res = take!(results)
        if !isnothing(res)
            write(w, join([res[x] for x = header], "\t"))
            write(w, "\n")
            flush(w)
        end
        next!(p)
        n = n - 1
    end

    close(w)
end