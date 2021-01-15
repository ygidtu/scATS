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

    export  BEDRecord, new_bed, ExtractedData, get_record_from_bam, run, setLevel

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
        # chrom = lines[2]
        # if !isnothing(chromosomes) && length(chromosomes) > 0
        #     if !haskey(chromosomes, chrom)
        #         if startswith(chrom, "chr")
        #             chrom = replace(chrom, "chr"=>"")
        #         else
        #             chrom = string("chr", chrom)
        #         end

        #         if !haskey(chromosomes, chrom)
        #             return BEDRecord("", 0, 0, "", "", "")
        #         end
        #     end
        # end

        # return BEDRecord(
        #     chrom, parse(Int64, lines[3]), parse(Int64, lines[4]),
        #     lines[6], lines[11], lines[1]
        # )
        
        if length(lines) >= 6
            return BEDRecord(
                replace(lines[1], r"chr" => ""), parse(Int64, lines[2]), parse(Int64, lines[3]),
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

    @with_kw struct ExtractedData
        utr::BEDRecord
        st_arr::Vector
        en_arr::Vector
    end

    function get_hash(data::ExtractedData)::String
        return string(
            join(sort(data.st_arr), ","), "|", 
            join(sort(data.en_arr), ","), "|",
            data.utr.end_pos - data.utr.start_pos
        )
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
        
        for record in eachoverlap(reader, utr.chrom, utr.start_pos:utr.end_pos)
            if !filter(record)
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

            # R2 end site, for Rn
            end_site = utr.strand == "-" ? BAM.leftposition(record) : BAM.rightposition(record)
            
            # calcualte R1 start site  based on SE mode and strand
            start_site = NaN
            if single_end_mode
                push!(st_arr, start_site)
            else
                start_site = BAM.nextposition(record)
                if utr.strand == "-"
                    start_site = start_site + BAM.alignlength(record)
                end
                push!(st_arr, start_site - utr_site)
            end
            push!(en_arr, end_site - utr_site)
        end

        if length(st_arr) > min_reads
            push!(
                res, 
                ExtractedData(utr, st_arr, en_arr)
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
        # inference with fixed parameters
        fixed_inference_flag::Bool = false,

        #single end mode
        single_end_mode::Bool = false,
        verbose::Bool = false,
        using_R::Bool = false,
        error_log=nothing
    )::Vector{String}
        res = Vector{String}()

        utr = data.utr
        en_arr = data.en_arr
        st_arr  = data.st_arr

        runner = using_R ? ATSMIX.fit_by_R : ATSMIX.fit

        temp_res = runner(
            n_max_ats, n_min_ats, 
            st_arr , en_arr; 
            L = utr.end_pos - utr.start_pos, 
            mu_f=mu_f, min_ws = min_ws, 
            sigma_f = sigma_f, max_beta=max_beta,
            fixed_inference_flag = fixed_inference_flag, 
            single_end_mode = single_end_mode, 
            verbose = verbose, error_log=error_log
        )
        if isnothing(temp_res.bic) || isinf(temp_res.bic)
            return res
        end

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

        return res
    end
end
