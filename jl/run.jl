
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
            help="max number of reads to test, to reduce memory usage, -1 mean to test all reads"
            arg_type = Int64
            default = -1
        "--seed"
            help="seed used to select reads to test, if there are too much reads"
            arg_type = Int64
            default = 42
        "--min-reads"
            help = "minimum reads to construct ATS"
            arg_type = Int
            default = 10
        # "--min-distance"
        #     help="minimum distance between to sites to test"
        #     arg_type = Int
        #     default = 10
        # "--max-distance"
        #     help="maximum distance between to sites to test"
        #     arg_type = Int
        #     default = 1000
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
    include(joinpath(@__DIR__, "src", "Extract.jl"))
end


function main(
    input_file::String, bam::String, output::String, distance::Int;
    mu_f::Int64=300, min_ws::Float64=0.01, min_pa_gap::Int64=100,
    max_beta::Int64=70, theta_step::Int64=9, sigma_f::Int64=50,
    using_R::Bool=false, verbose::Bool=false,
    max_reads_to_test::Int64=-1, seed::Int=42,
    process::Int=1, min_reads::Int=10, 
    # min_distance::Int=10, max_distance::Int=500
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

    reader = open(BAM.Reader, bam)
    chromosomes = Dict{String, String}()
    for h = findall(header(reader), "SQ")
        chromosomes[h["SN"]] = h["LN"]
    end
    close(reader)

    # println(chromosomes)

    beds = Dict{String, Vector}()
    r = open(input_file, "r")
    seekend(r)
    fileSize = position(r)

    seekstart(r)
    p = Progress(fileSize, 1)   # minimum update interval: 1 second
    while !eof(r)
        try
            temp_bed = Extract.new_bed(readline(r), chromosomes=chromosomes)

            if temp_bed.chrom == ""
                continue
            end

            temp = get(beds, temp_bed.score, Vector())

            if length(temp) == 0 || temp[length(temp)].name != temp_bed.name
                push!(temp, temp_bed)
            end
            beds[temp_bed.score] = temp
            # if length(beds) > 10
            #     break
            # end
        catch e
            warn(logger, string("error while reading promoters: ", e))
            continue
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
    res = @showprogress "Reading..." pmap(collect(values(beds))) do b
        temp = Extract.get_record_from_bam(
            bam, b, 
            distance=distance, 
            min_reads=min_reads, 
            # min_distance=min_distance,
            # max_distance=max_distance,
            max_number_of_sites=max_reads_to_test,
            seed=seed
        )
        return temp
    end

    # extract reads and drop duplicates
    temp_data = Dict{String, Vector{Extract.ExtractedData}}()
    total = 0
    for t = res
        for r = t
            total += 1
            k = Extract.get_hash(r)

            temp_vec = get(temp_data, k, Vector())
            push!(temp_vec, r)
            temp_data[k] = temp_vec
        end
    end

    info(logger, string("compressed: ", length(temp_data), "; total: ", total))
    # calculate
    res = @showprogress "Computing..." pmap(collect(values(temp_data))) do data
        return Extract.run(data)
    end

    open(output, "w+") do w
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
    args["input"], args["bam"], args["output"], args["distance"],
    mu_f = args["mu-f"], sigma_f = args["sigma-f"],
    min_ws = args["min-ws"], min_pa_gap=args["min-pa-gap"],
    max_beta=args["max-beta"], theta_step=args["theta-step"],
    using_R=args["using-R"], verbose=args["verbose"],
    max_reads_to_test=args["max-reads-to-test"], seed=args["seed"],
    process=args["process"], min_reads=args["min-reads"], 
    # min_distance=args["min-distance"], max_distance=args["max-distance"]
)