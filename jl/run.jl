
using ArgParse
using Distributed


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
            default = -500
        "--using-R"
            help="whether to use R ksmooth"
            action = :store_true
        "--process", "-p"
            help = "How many processes to use"
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
        "--min-pa-gap"
            help = "min pa gap"
            arg_type = Int64
            default = 100
        "--max-beta"
            help = "maximum beta"
            arg_type = Int64
            default = 70
        "--theta-step"
            help="theta step"
            arg_type = Int64
            default = 9
        "--max-reads-to-test"
            help="max number of reads to test, to reduce memory usage"
            arg_type = Int64
            default = 200
        "--seed"
            help="seed used to select reads to test, if there are too much reads"
            arg_type = Int64
            default = 42
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
    # include(joinpath(@__DIR__, "src", "APA.jl"))
    include(joinpath(@__DIR__, "src", "Extract.jl"))
    #=
    struct Data
        BC::String
        Chrom::String
        UTRStart::Int64
        UTREnd::Int64
        ReadStart::Int64
        ReadEnd::Int64
        ReadLabel::String
        StartLocInUTR::Union{Int64, Float64}
        LenInUTR::Union{Int64, Float64}
        LenInPA::Union{Int64, Float64}
        LenPA::Union{Int64, Float64}
        PASite::Union{Int64, Float64}
    end

    function newData(data::Dict{String, String})::Data
        start_loc_in_utr = get(data, "StartLocInUTR", nothing)
        if !isnothing(start_loc_in_utr)
            start_loc_in_utr = parse(Int64, start_loc_in_utr)
        end
        len_in_utr = get(data, "LenInUTR", "NA")
        if len_in_utr != "NA"
            len_in_utr = parse(Int64, len_in_utr)
        else
            len_in_utr = NaN
        end
        len_in_pa = get(data, "LenInPA", "NA")
        if len_in_pa != "NA"
            len_in_pa = parse(Int64, len_in_pa)
        else
            len_in_pa = NaN
        end
        len_pa = get(data, "LenPA", "NA")
        if len_pa != "NA"
            len_pa = parse(Int64, len_pa)
        else
            len_pa = NaN
        end
        pa_site = get(data, "PASite", "NA")
        if pa_site != "NA"
            pa_site = parse(Int64, pa_site)
        else
            pa_site = NaN
        end

        return Data(
            data["BC"], data["Chrom"], 
            parse(Int64, data["UTRStart"]), parse(Int64, data["UTREnd"]),
            parse(Int64, data["ReadStart"]), parse(Int64, data["ReadEnd"]),
            data["ReadLabel"], start_loc_in_utr, len_in_utr,
            len_in_pa, len_pa, pa_site
        )
  
    end

    function runAPA(
        data::Vector;
        mu_f::Int64=300, min_ws::Float64=0.01, min_pa_gap::Int64=100,
        max_beta::Int64=70, theta_step::Int64=9, sigma_f::Int64=50,
        using_R::Bool=false, verbose::Bool=false
    )::Vector{String}
        data = [newData(i) for i = data]
        res = Vector{String}()
        for i = data
            #=
            r1_utr_st_arr=df.StartLocInUTR # start location of each read on UTR part

            r1_utr_st_arr = r1_utr_st_arr
            r1_len_arr=df.LenInUTR # length of each read on UTR part
            polya_len_arr=[NaN for _ in df.LenPA] # length of polyA
            # pa_site_arr=df$PASite, # pa site locations of reads
            L = mean(df.UTREnd - df.UTRStart)

            r1_utr_st_arr = r1_utr_st_arr .+ L
            =#
            temp_data = data[findall(x -> x.UTRStart == i.UTRStart && x.UTREnd == i.UTREnd && x.BC == i.BC, data)]
            if length(temp_data) < 2
                continue
            end
            r1_utr_st_arr, r1_len_arr, polya_len_arr, Ls = Vector{Int64}(), Vector{Int64}(), Vector{Union{Int64, Float64}}(), Vector{Int64}()
            reads_start, reads_end = Vector{String}(), Vector{String}()
            for j = temp_data
                push!(r1_utr_st_arr, j.StartLocInUTR)
                push!(r1_len_arr, j.LenInUTR)
                push!(polya_len_arr, j.LenPA)
                push!(Ls, j.UTREnd - j.UTRStart)
                push!(reads_start, string(j.ReadStart))
                push!(reads_end, string(j.ReadEnd))
            end
            L = mean(Ls)

            try
                self = APA.new(
                    5, 1,
                    r1_utr_st_arr .+ L,
                    r1_len_arr,
                    polya_len_arr,
                    polya_len_arr,
                    polya_len_arr,
                    trunc(Int64, maximum(r1_len_arr) + maximum(r1_utr_st_arr) + 300),
                    mu_f = mu_f,
                    sigma_f = sigma_f,
                    min_ws = min_ws,
                    min_pa_gap=min_pa_gap,
                    max_beta=max_beta, 
                    theta_step=theta_step,
                    verbose = verbose
                )

                temp_res = APA.fit(self, using_R=using_R)

                push!(res,
                    Formatting.format(
                        FormatExpr("{}\t{}:{}-{}\t{}\t{}\t{}"), 
                        i.BC, i.Chrom, i.UTRStart, i.UTREnd, 
                        join(reads_start, ","), join(reads_end, ","), 
                        string(temp_res)
                    )
                )
            catch e
                println(e)
                println(r1_utr_st_arr)
                println(r1_len_arr)
                println(L)
                continue
            end

        end
        return res
    end
    =#
end


function main(
    input_file::String, bam::String, output::String, distance::Int;
    mu_f::Int64=300, min_ws::Float64=0.01, min_pa_gap::Int64=100,
    max_beta::Int64=70, theta_step::Int64=9, sigma_f::Int64=50,
    using_R::Bool=false, verbose::Bool=false,
    max_reads_to_test::Int64=200, seed::Int=42
)
    logger = Memento.config!("debug"; fmt="[{date} - {level} | {name}]: {msg}")
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

    #=
    HEADER = ["BC", "Chrom", "UTRStart", "UTREnd", "ReadStart", "ReadEnd", "ReadLabel", "StartLocInUTR", "LenInUTR", "LenInPA", "LenPA", "PASite"]
    processed = ""

    temp_data = Vector()

    open(input_file,"r") do file
        seekend(file)
        fileSize = position(file)
        seekstart(file)
        p = Progress(fileSize, 1)   # minimum update interval: 1 second
        open(output, "w+") do w
            while !eof(file)
                line = split(strip(readline(file)), "\t")
                update!(p, position(file))
                data = Dict(i=>string(j) for (i, j) in zip(HEADER, line))
                # println(data)
                if parse(Int, get(data, "StartLocInUTR", -Inf)) < distance
                    continue
                end
                
                if  processe != data["Chrom"] && length(temp_data) > 100
                    # process previous data
                    if length(temp_data) > 1

                        res = pmap(collect(keys(temp_data))) do p
                            return runAPA(
                                temp_data[p],
                                mu_f = mu_f, sigma_f = sigma_f,
                                min_ws = min_ws, min_pa_gap = min_pa_gap,
                                max_beta = max_beta, theta_step = theta_step,
                                using_R = using_R, verbose = verbose
                            )
                        end

                        for i = res
                            write(w, string(join(i, "\n"), "\n"))
                        end
                        flush(w)
                    end

                    # reset collected data
                    processed = data["Chrom"]
                    temp_data = [data]
                else
                    # if haskey(temp_data, data["BC"])
                    #     push!(temp_data[data["BC"]], data)
                    # else
                    #     temp_data[data["BC"]] = [data]
                    # end
                    push!(temp_data, data)
                end
            end
            close(w)
        end

        close(file)
    end
    =#

    beds = Vector()
    r = open(input_file, "r")
    seekend(r)
    fileSize = position(r)

    seekstart(r)
    p = Progress(fileSize, 1)   # minimum update interval: 1 second
    while !eof(r)
        push!(beds, Extract.new_bed(readline(r)))
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
    # res = @showprogress pmap(beds) do b
    #     # println(b)
    #     return Extract.get_record_from_bam(
    #         bam, b, distance=distance,
    #         mu_f = mu_f, sigma_f = sigma_f,
    #         min_ws = min_ws, min_pa_gap = min_pa_gap,
    #         max_beta = max_beta, theta_step = theta_step,
    #         using_R = using_R, verbose = verbose
    #     )
    #     # println(temp)
    #     # return temp
    # end

    open(output, "w+") do w
        chunk = 1000
        pbar = Progress(length(beds), 1)
        for i in 1:chunk:length(beds)

            # res = Vector()
            # for b = beds[i:min(i+chunk, length(beds))]
            res = pmap(beds[i:min(i+chunk, length(beds))]) do b
                temp = Extract.get_record_from_bam(
                    bam, b, distance=distance,
                    mu_f = mu_f, sigma_f = sigma_f,
                    min_ws = min_ws, min_pa_gap = min_pa_gap,
                    max_beta = max_beta, theta_step = theta_step,
                    using_R = using_R, verbose = verbose,
                    max_number_of_sites=max_reads_to_test, seed=seed
                )
                # push!(res, temp)
                return temp
            end


            for r = res
                if r != ""
                    write(w, string(r, "\n"))
                end
            end
        end
        close(w)
    end
end


main(
    args["input"], args["bam"], args["output"], args["distance"],
    mu_f = args["mu-f"], sigma_f = args["sigma-f"],
    min_ws = args["min-ws"], min_pa_gap=args["min-pa-gap"],
    max_beta=args["max-beta"], theta_step=args["theta-step"],
    using_R=args["using-R"], verbose=args["verbose"],
    max_reads_to_test=args["max-reads-to-test"], seed=args["seed"]
)