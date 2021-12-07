
include(joinpath(@__DIR__, "src", "bam.jl"))
include(joinpath(@__DIR__, "src", "core.jl"))
include(joinpath(@__DIR__, "src", "genomic.jl"))


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

    l = ReentrantLock()

    # Counting by multi-threads
    p = Progress(length(utrs), 1, "Computing...")

    w = open(output, "w+")

    # write header
    header = ATSModel.getHeader()
    write(w, string(join(header, "\t"), "\n"))
    flush(w)

    Threads.@threads for utr = utrs
        st_arr = Bam.load_reads(bams, utr)

        if length(st_arr) < args["min-reads"]
            continue
        end

        param = ATSModel.init(
            st_arr, 
            utr,
            n_max_ats=args["n-max-ats"], 
            n_min_ats=args["n-min-ats"], 
            min_ws=args["min-ws"], 
            max_unif_ws=args["max-unif-ws"], 
            fixed_inference_flag=args["fixed-inference"],
            step_size=args["step-size"],
            nround=args["nround"],
            max_beta=args["max-beta"]
        )

        res = ATSModel.run(param)
        res = ATSModel.toDict(res)

        lock(l)
        write(w, join([res[x] for x = header], "\t"))
        write(w, "\n")
        flush(w)
        unlock(l)

        next!(p)
    end

    close(w)
end