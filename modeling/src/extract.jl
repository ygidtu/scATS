module Extract
    using BioAlignments
    using Compat: @__MODULE__  # requires a minimum of Compat 0.26. Not required on Julia 0.7
    using FilePathsBase
    using Formatting
    using GenomicFeatures
    using Memento
    using Parameters
    using ProgressMeter
    using Random
    using StatsBase
    using XAM

    include(joinpath(@__DIR__, "ATS.jl"))

    # Create our module level LOGGER (this will get precompiled)
    const LOGGER = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg} | {stacktrace}")
    # const LOGGER = getlogger(@__MODULE__)

    # Register the module level LOGGER at runtime so that folks can access the LOGGER via `get_LOGGER(MyModule)`
    # NOTE: If this line is not included then the precompiled `MyModule.LOGGER` won't be registered at runtime.
    function __init__()
        Memento.register(LOGGER)
    end

    export  BEDRecord, new_bed, ExtractedData, get_record_from_bam, run, setLevel, get_bed_short

    function setLevel(level::String)
        setlevel!(LOGGER, level)
        ATSMIX.setLevel(level)
    end

    @with_kw struct BEDRecord
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
    
    function new_bed(line::String)::BEDRecord  # ; chromosomes::Dict{String, String}=nothing
        lines = split(strip(line), "\t")
        
        if length(lines) >= 6
            return BEDRecord(
                lines[1], parse(Int64, lines[2]), parse(Int64, lines[3]),
                lines[6], lines[4], strip(lines[5])
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

    function get_bed_short(bed::BEDRecord)::String
        return Formatting.format(
            FormatExpr("{}:{}-{}:{}"),
            self.chrom, self.start_pos, self.end_pos, self.strand
        )
    end

    @with_kw struct ExtractedData
        utr::BEDRecord
        st_arr::Vector
        en_arr::Vector
        real_st::Vector
        real_en::Vector
    end

    Base.show(io::IO, self::ExtractedData) = print(
        io,
        Formatting.format(
            FormatExpr("utr: {}\nst_arr: {}\nen_arr: {}\n"), # real_st: {}\nreal_en: {}
            string(self.utr),
            join(map(string, self.st_arr), ","),
            join(map(string, self.en_arr), ","),
            # join(map(string, self.real_st), ","),
            # join(map(string, self.real_en), ",")
        )
    )

    function get_utr(data::ExtractedData)::String
        return string(data.utr.chrom, "_", data.utr.start_pos, "_", data.utr.end_pos, "_", data.utr.strand)
    end

    function filter(record)::Bool
        #=
        Filter low quality reads
        =#

        if !BAM.ismapped(record)
            return false
        end

        auxdata = Dict(BAM.auxdata(record))
        if haskey(auxdata, "NH") && auxdata["NH"] > 1
            return false
        end

        if BAM.flag(record) & SAM.FLAG_QCFAIL > 0
            return false
        end

        return true # haskey(auxdata, "CB") && haskey(auxdata, "UB")
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
        path::String, 
        utr::BEDRecord;
        single_end_mode::Bool = false,
        min_reads::Int=0,
        distance::Int=1500
    )::Vector{ExtractedData}
        bai = string(path, ".bai")
        promoter_sites = Dict{BEDRecord, Vector}()
        reader = open(BAM.Reader, path, index=bai)

        sites = Dict{}
        res = Vector()

        st_arr = Vector()
        en_arr = Vector()
        real_st, real_en = Vector(), Vector()
        
        for record in eachoverlap(reader, utr.chrom, utr.start_pos:utr.end_pos)
            if !filter(record)
                continue
            end

            if BAM.flag(record) & SAM.FLAG_PROPER_PAIR == 0
                continue
            end                
            
            # Only kept R2
            if BAM.flag(record) & SAM.FLAG_READ1 > 0
                continue
            end

            if !single_end_mode && !BAM.isnextmapped(record)
                continue
            end

            utr_site = utr.strand == "-" ? utr.start_pos : utr.end_pos

            # R2 needs locate in UTR
            if utr.start_pos > BAM.leftposition(record) || utr.end_pos < BAM.rightposition(record)
                continue
            end

            # R1 needs locate in UTR
            if !single_end_mode && (utr.start_pos > BAM.nextposition(record) || utr.end_pos < BAM.nextposition(record) + BAM.alignlength(record))
                continue
            end

            # read strand
            strand = determine_strand(record)

            if strand != utr.strand
                continue
            end

            # R2 end site, for Rn
            end_site = strand == "-" ? BAM.leftposition(record) : BAM.rightposition(record)

            # calcualte R1 start site  based on SE mode and strand
            start_site = NaN
            if single_end_mode
                push!(st_arr, start_site)
            else
                start_site = BAM.nextposition(record)
                if strand == "-"
                    start_site = start_site + BAM.alignlength(record)
                end

                push!(st_arr, start_site - utr_site)
                push!(real_st, start_site)
            end
            push!(en_arr, end_site - utr_site)
            push!(real_en, end_site)
        end

        if length(st_arr) > min_reads
            push!(
                res, 
                ExtractedData(utr, st_arr, en_arr, real_st, real_en)
            )
        end
        close(reader)

        return res
    end

    function run(
        data::ExtractedData;
        # maximum number of ATS sites
        n_max_ats::Int=5, 
        # minimum number of ATS sites
        n_min_ats::Int=1, 
        # fragment size information
        # fragment length mean
        mu_f::Int = 350, 
        # fragment length standard deviation
        sigma_f::Int = 50, 

        # pa site information
        # minimum weight of ATS site
        min_ws::Float64 = 0.01, 
        # maximum std for ATS site
        max_beta::Float64 = 50.0,
        seed::Int=2,
        # inference with fixed parameters
        fixed_inference_flag::Bool = false,

        #single end mode
        single_end_mode::Bool = false,
        verbose::Bool = false,
        using_R::Bool = false,
        error_log=nothing,
        debug::Bool = false,
        density=nothing
    )::Vector
        res = Vector()

        utr = data.utr
        en_arr = data.en_arr
        st_arr  = data.st_arr

        runner = using_R ? ATSMIX.fit_by_R : ATSMIX.fit

        temp_res = runner(
            n_max_ats, n_min_ats, 
            st_arr , en_arr; 
            L = abs(utr.end_pos - utr.start_pos), 
            mu_f=mu_f, min_ws = min_ws, 
            sigma_f = sigma_f, max_beta=max_beta,
            fixed_inference_flag = fixed_inference_flag, 
            single_end_mode = single_end_mode, 
            verbose = verbose, error_log=error_log,
            seed=seed, debug=density
        )
        if isnothing(temp_res.bic) || isinf(temp_res.bic)
            return res
        end

        if debug
            push!(res, temp_res)
        else
            push!(
                res, 
                Formatting.format(
                    FormatExpr("{}:{}-{}:{}\t{}\t{}\t{}"), 
                    data.utr.chrom, data.utr.start_pos, data.utr.end_pos, data.utr.strand,
                    join(map(string, data.st_arr), ","),
                    join(map(string, data.en_arr), ","),
                    string(temp_res)
                )
            )
        end

        
        return res
    end
end
