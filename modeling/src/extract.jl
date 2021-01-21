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
    const LOGGER = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")
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
                replace(lines[1], "chr"=>""), parse(Int64, lines[2]), parse(Int64, lines[3]),
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
            FormatExpr("utr: {}\nst_arr: {}\nen_arr: {}"), # \nreal_st: {}\nreal_en: {}
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


    function read_from_bam(path::String,  utr::BEDRecord; single_end_mode::Bool = false, min_reads::Int=0)::Vector{ExtractedData}
        bai = string(path, ".bai")
        reader = open(BAM.Reader, path, index=bai)
        sites = Dict{}
        res = Vector()

        st_arr, en_arr = Vector(), Vector()
        real_st, real_en = Vector(), Vector()
        r1, r2 = Dict(), Dict()
        
        for record in eachoverlap(reader, utr.chrom, utr.start_pos:utr.end_pos)
            if !filter(record)
                continue
            end

            if BAM.flag(record) & SAM.FLAG_PROPER_PAIR == 0
                continue
            end
            
            if !single_end_mode && !BAM.isnextmapped(record)
                continue
            end
            
            # Only kept R2
            if BAM.flag(record) & SAM.FLAG_READ1 > 0
                r1[BAM.tempname(record)] = Dict(
                    "start"=>BAM.leftposition(record), 
                    "end"=>BAM.rightposition(record)
                )
                continue
            end

            # R2 needs locate in UTR
            if utr.start_pos > BAM.leftposition(record) || utr.end_pos < BAM.rightposition(record)
                continue
            end

            # read strand
            strand = determine_strand(record)

            if strand != utr.strand
                continue
            end

            r2[BAM.tempname(record)] = Dict(
                "start"=>BAM.leftposition(record), 
                "end"=>BAM.rightposition(record)
            )           
        end

        # now filter and check r1
        utr_site = utr.strand == "-" ? utr.start_pos : utr.end_pos
        for (name, site) = r2

            r1_pos = get(r1, name, nothing)

            if isnothing(r1_pos)
                continue
            end

            # R1 needs locate in UTR
            if !single_end_mode && (utr.start_pos > r1_pos["start"] || utr.end_pos < r1_pos["end"])
                continue
            end

            # R2 end site, for Rn
            end_site = strand == "-" ? site["start"] : site["end"]

            push!(en_arr, end_site - utr_site)
            push!(real_en, end_site)

            # calcualte R1 start site  based on SE mode and strand
            start_site = NaN
            if single_end_mode
                push!(st_arr, start_site)
            else
                start_site = utr.strand == "+" ? r1_pos["start"] : r1_pos["end"]
                # start_site = r1_pos["start"]

                push!(st_arr, start_site - utr_site)
                push!(real_st, start_site)
            end
        end

        if length(st_arr) > min_reads
            push!(
                res, 
                ExtractedData(utr, st_arr, en_arr, real_st, real_en) # 
            )
        end
        close(reader)

        return res
    end

    function read_from_bam_cage(path::String,  utr::BEDRecord; min_reads::Int=0)
        bai = string(path, ".bai")
        reader = open(BAM.Reader, path, index=bai)
        sites = Dict{}
        res = Vector()

        st_arr, en_arr = Vector(), Vector()
        real_st, real_en = Vector(), Vector()
        utr_site = utr.strand == "-" ? utr.start_pos : utr.end_pos
        for record in eachoverlap(reader, utr.chrom, utr.start_pos:utr.end_pos)
            if !filter(record)
                continue
            end

            # Only kept R1
            if BAM.flag(record) & SAM.FLAG_READ2 > 0
                continue
            end

            if utr.start_pos > BAM.leftposition(record) || utr.end_pos < BAM.rightposition(record)
                continue
            end

            # read strand
            strand = determine_strand(record)

            if strand != utr.strand
                continue
            end

            start_site = utr.strand == "+" ? BAM.leftposition(record) : BAM.rightposition(record)
            push!(st_arr, start_site - utr_site)
            push!(real_st, start_site)
        end

        if length(st_arr) > min_reads
            push!(
                res, 
                ExtractedData(utr, st_arr, en_arr, real_st, real_en) # 
            )
        end
        close(reader)
        return res
    end


    function get_record_from_bam(
        path::String, 
        utr::BEDRecord;
        single_end_mode::Bool = false,
        min_reads::Int=0,
        cage_mode::Bool = false
    )::Vector{ExtractedData}

        if cage_mode
            return read_from_bam_cage(path, utr, min_reads=min_reads)
        end
        return read_from_bam(path, utr, single_end_mode=single_end_mode, min_reads=min_reads)
        
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
        cage_mode::Bool = false,
        using_R::Bool = false,
        error_log=nothing,
        debug::Bool = false,
    )::Vector
        res = Vector()

        utr = data.utr
        en_arr = data.en_arr
        st_arr  = data.st_arr

        temp_res = ATSMIX.fit(
            n_max_ats, n_min_ats, 
            st_arr , en_arr; 
            L = abs(utr.end_pos - utr.start_pos), 
            mu_f=mu_f, min_ws = min_ws, 
            sigma_f = sigma_f, max_beta=max_beta,
            fixed_inference_flag = fixed_inference_flag, 
            single_end_mode = single_end_mode, 
            cage_mode=cage_mode,
            error_log=error_log,
            seed=seed, using_R = using_R
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
                    # join(map(string, data.real_st), ","),
                    # join(map(string, data.real_en), ","),
                    string(temp_res)
                )
            )
        end

        return res
    end
end
