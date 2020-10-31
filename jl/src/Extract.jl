module Extract
    using BioAlignments
    using Compat: @__MODULE__  # requires a minimum of Compat 0.26. Not required on Julia 0.7
    using FilePathsBase
    using Formatting
    using GenomicFeatures
    using Memento
    using ProgressMeter
    using Random
    using StatsBase

    include(joinpath(@__DIR__, "APA.jl"))

    # Create our module level LOGGER (this will get precompiled)
    Memento.config!("debug"; fmt="[{date} - {level} | {name}]: {msg}")
    const LOGGER = getlogger(@__MODULE__)

    # Register the module level LOGGER at runtime so that folks can access the LOGGER via `get_LOGGER(MyModule)`
    # NOTE: If this line is not included then the precompiled `MyModule.LOGGER` won't be registered at runtime.
    function __init__()
        Memento.register(LOGGER)
    end

    struct BEDRecord
        chrom::String
        start_pos::Int64
        end_pos::Int64
        strand::String
        score::String
        name::String 
    end

    Base.show(io::IO, self::BEDRecord) = print(
        io,
        Formatting.format(
            FormatExpr("{}:{}-{}:{}\t{}\t{}"),
            self.chrom, self.start_pos, self.end_pos, self.strand,
            self.name, self.score
        )
    )
    
    function new_bed(line::String)::BEDRecord
        lines = split(strip(line), "\t")
        
        if length(lines) >= 6
            return BEDRecord(
                lines[1], parse(Int64, lines[2]), parse(Int64, lines[3]),
                lines[6], lines[4], lines[5]
            )
        elseif length(lines) == 3
            return BEDRecord(
                lines[1], parse(Int64, lines[2]), parse(Int64, lines[3]),
                "", "", ""
            )
        elseif length(lines) == 4
            return BEDRecord(
                lines[1], parse(Int64, lines[2]), parse(Int64, lines[3]),
                strip(lines[4]), "", ""
            )
        else
            error(LOGGER, string("the number of columns mismatch: ", lines))
            exit(1)
        end
    end

    function filter(record)::Bool
        #=
        Filter low quality reads
        =#

        auxdata = Dict(BAM.auxdata(record))
        if haskey(auxdata, "NH") && auxdata["NH"] > 1
            return false
        end

        if BAM.flag(record) & SAM.FLAG_DUP > 0 || BAM.flag(record) & SAM.FLAG_QCFAIL > 0
            return false
        end

        return haskey(auxdata, "CB") && haskey(auxdata, "UB")
    end

    function determine_strand(record)::String
        strand = "*"
        flag = BAM.flag(record)
        if flag & SAM.FLAG_READ1 > 0
            if flag & SAM.FLAG_REVERSE > 0
                strand = "-"
            else
                strand = "+"
            end
        elseif flag & SAM.FLAG_READ2 > 0
            if flag & SAM.FLAG_REVERSE > 0
                strand = "+"
            else
                strand = "-"
            end
        end
        return strand
    end

    function get_record_from_bam(
        path::String, utr::BEDRecord;
        max_number_of_sites::Int64=200, seed::Int64=42,
        distance::Int64=-500, mu_f::Int64=300, min_ws::Float64=0.01, 
        min_pa_gap::Int64=100, max_beta::Int64=70, theta_step::Int64=9, 
        sigma_f::Int64=50, using_R::Bool=false, verbose::Bool=true
    )::String
        data, polya_len_arr = Vector{Vector{Int64}}(), Vector{Union{Int64, Float64}}()
        reads_start, reads_end = Vector{String}(), Vector{String}()

        bai = string(path, ".bai")
        # println(string(utr, " start"))
        
        reader = open(BAM.Reader, path, index=bai)
        try
            for record in eachoverlap(reader, utr.chrom, utr.start_pos:utr.end_pos)
                if !filter(record)
                    continue
                end

                strand = utr.strand

                utr_site = utr.end_pos
                start_site = BAM.position(record)
                if strand != "+"
                    utr_site = utr.start_pos
                    start_site = BAM.rightposition(record)
                end

                if BAM.position(record) <= utr_site <= BAM.rightposition(record) && start_site - utr_site > distance
                    t1 = start_site-utr_site
                    t2 = min(BAM.rightposition(record), utr.end_pos) - max(BAM.position(record), utr.start_pos)

                    push!(data, [t1, t2])
                    push!(polya_len_arr, NaN)
                    push!(reads_start, string(BAM.position(record)))
                    push!(reads_end, string(BAM.rightposition(record)))
                end
            end
            close(reader)
        catch e
            close(reader)
            warn(LOGGER, string(utr))
            error(LOGGER, e)
            return ""
        end

        # data = unique(data)
        # maxmium used 200 (default) sites to save memory
        if length(data) > max_number_of_sites
            Random.seed!(seed)
            data = Random.shuffle(data)
            data = data[1:max_number_of_sites]
            polya_len_arr = polya_len_arr[1:max_number_of_sites]
        end

        if length(data) > 1
            try
                r1_utr_st_arr = [x[1] for x = data]
                r1_len_arr = [x[2] for x = data]
                self = APA.new(
                    5, 1,
                    r1_utr_st_arr .+ (utr.end_pos - utr.start_pos),
                    r1_len_arr,
                    polya_len_arr[1:length(data)],
                    polya_len_arr[1:length(data)],
                    polya_len_arr[1:length(data)],
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

                return Formatting.format(
                    FormatExpr("{}:{}-{}\t{}\t{}\t{}"), 
                    utr.chrom, utr.start_pos, utr.end_pos, 
                    join(reads_start, ","), join(reads_end, ","), 
                    string(temp_res)
                )
            catch e
                warn(LOGGER, e)
                if verbose
                    debug(LOGGER, join([string(x[1]) for x = data], ","))
                    debug(LOGGER, join([string(x[2]) for x = data], ","))
                    debug(LOGGER, string(utr.end_pos - utr.start_pos))
                end
            end
        end

        return ""
    end
end