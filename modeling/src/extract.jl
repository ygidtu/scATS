module Extract
    using BGZFStreams
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
    include(joinpath(@__DIR__, "genomic.jl"))

    # Create our module level LOGGER (this will get precompiled)
    const LOGGER = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")
    # const LOGGER = getlogger(@__MODULE__)

    # Register the module level LOGGER at runtime so that folks can access the LOGGER via `get_LOGGER(MyModule)`
    # NOTE: If this line is not included then the precompiled `MyModule.LOGGER` won't be registered at runtime.
    function __init__()
        Memento.register(LOGGER)
    end

    export  TestData, get_bed_short, get_record_from_bam, get_record_from_bam_transcript, run, setLevel

    function setLevel(level::String)
        setlevel!(LOGGER, level)
        ATSMIX.setLevel(level)
    end

    @with_kw struct TestData
        utr
        st_arr::Vector
        en_arr::Vector
        real_st::Vector = []
        real_en::Vector = []
        exon_coord::Dict{Int, Int} = Dict()
    end

    Base.show(io::IO, self::TestData) = print(
        io,
        Formatting.format(
            FormatExpr("utr: {}\nst_arr: {}\nen_arr: {}"), # \nreal_st: {}\nreal_en: {}
            string(self.utr),
            join(map(string, self.st_arr), ","),
            join(map(string, self.en_arr), ","),
        )
    )

    function get_utr(data::TestData)::String
        return string(data.utr.chrom, "_", data.utr.start_pos, "_", data.utr.end_pos, "_", data.utr.strand)
    end

    @with_kw struct Reads
        start_pos::Int
        end_pos::Int
        intron_site::Vector{Int}
        source::Union{String, Nothing} = nothing
    end

    @with_kw struct ExtractData
        utr
        R1::Vector{Reads}
        R2::Vector{Reads}
    end

    function toDict(self)::Dict
        if typeof(self) == ExtractData
            return Dict(
                toDict(self.utr) => Dict(
                    Symbol("R1") => toDict(self.R1),
                    Symbol("R2") => toDict(self.R2)
                )
            )
        end

        return Dict(x => getfield(self, x) for x = fieldnames(typeof(self)))
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


    function cigar_to_intron(record)::Vector
        res = Vector()
        cigar = BAM.cigar(record)
        init = BAM.leftposition(record)
        dist = ""
        skip_code = Dict('I' => 0, 'D' => 0, 'H' => 0, 'S' => 0)
        for i = cigar
            try
                parse(Int, i)
                dist = string(dist, i)             
            catch
                if !haskey(skip_code, i)
                    init += parse(Int, dist)
                end
            
                if i == 'N'
                    push!(res, init - parse(Int, dist))
                    push!(res, init)
                end
                dist = ""
            end
        end
        return res
    end

    function read_from_bam(path::Vector,  utr; single_end_mode::Bool = false, min_reads::Int=0)::Union{TestData, ExtractData, Nothing}
        st_arr, en_arr = Vector(), Vector()
        real_st, real_en = Vector(), Vector()
        r1, r2 = Dict(), Dict()

        for p = path
            key = basename(Path(p))
            bai = string(p, ".bai")

            reader = open(BAM.Reader, p, index=bai)

            for record in eachoverlap(reader, utr.Chrom, utr.Start:utr.End)
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
                    r1[string(key, "|", BAM.tempname(record))] = record
                    continue
                end

                # R2 needs locate in UTR
                if utr.Start > BAM.leftposition(record) || utr.End < BAM.rightposition(record)
                    continue
                end

                # read strand
                strand = determine_strand(record)
                if strand != utr.Strand
                    continue
                end

                r2[string(key, "|", BAM.tempname(record))] = record
            end
            close(reader)
        end

        utr_site = utr.Strand == "+" ? utr.Start : utr.End
        for (name, site) = r2
            r1_pos = get(r1, name, nothing)

            if isnothing(r1_pos)
                continue
            end

            # R1 needs locate in UTR
            if !single_end_mode && (utr.Start > BAM.leftposition(r1_pos) || utr.End < BAM.rightposition(r1_pos))
                continue
            end

            # R2 end site, for Rn
            if is_extract
                push!(st_arr, Reads(
                    start_pos = BAM.leftposition(r1_pos),
                    end_pos = BAM.rightposition(r1_pos),
                    intron_site = cigar_to_intron(r1_pos),
                    source = split(name, "|")[1]
                ))

                push!(en_arr, Reads(
                    start_pos = BAM.leftposition(site),
                    end_pos = BAM.rightposition(site),
                    intron_site = cigar_to_intron(site),
                    source = split(name, "|")[1]
                ))
            else
                end_site = strand == "-" ? BAM.leftposition(site) : BAM.rightposition(site)

                push!(en_arr, end_site - utr_site)
                push!(real_en, end_site)

                # calcualte R1 start site  based on SE mode and strand
                start_site = NaN
                if single_end_mode
                    push!(st_arr, start_site)
                else
                    start_site = utr.Strand == "+" ? BAM.leftposition(r1_pos) : BAM.rightposition(r1_pos)
                    
                    push!(st_arr, start_site - utr_site)
                    push!(real_st, start_site)
                end
            end
        end

        if is_extract
            return ExtractData(utr = utr, R1 = st_arr, R2 = en_arr)
        elseif length(st_arr) > min_reads
            return TestData(utr = utr, st_arr = st_arr, en_arr = en_arr, real_st = real_st, real_en = real_en)
        end

        return nothing
    end

    function read_from_bam_cage(path::Vector, utr; min_reads::Int=0)::Union{TestData, Nothing}

        st_arr, en_arr = Vector(), Vector()
        real_st, real_en = Vector(), Vector()
        r1, r2 = Dict(), Dict()

        for p = path
            key = basename(Path(p))
            bai = string(p, ".bai")
 
            reader = open(BAM.Reader, p, index=bai)
            sites = Dict{}

            utr_site = utr.Strand == "+" ? utr.Start : utr.End
            for record in eachoverlap(reader, utr.Chrom, utr.Start:utr.End)
                if !filter(record)
                    continue
                end

                # Only kept R1
                if BAM.flag(record) & SAM.FLAG_READ2 > 0
                    r2[string(key, "|", BAM.tempname(record))] = Dict(
                        "start"=>BAM.leftposition(record), 
                        "end"=>BAM.rightposition(record),
                    )
                    continue
                end

                if utr.Start > BAM.leftposition(record) || utr.End < BAM.rightposition(record)
                    continue
                end

                # read strand
                strand = determine_strand(record)

                if strand != utr.Strand
                    continue
                end

                r1[string(key, "|", BAM.tempname(record))] = Dict(
                    "start"=>BAM.leftposition(record), 
                    "end"=>BAM.rightposition(record),
                )
            end
            close(reader)
        end

        utr_site = utr.Strand == "+" ? utr.Start : utr.End
        for (name, r1_pos) = r1
            r2_pos = get(r2, name, nothing)

            if isnothing(r2_pos)
                continue
            end

            # R1 needs locate in UTR
            if !single_end_mode && (utr.Start > r1_pos["start"] || utr.End < r1_pos["end"])
                continue
            end

            # R2 end site, for Rn
            end_site = strand == "-" ? r2_pos["start"] : r2_pos["end"]

            push!(en_arr, end_site - utr_site)
            push!(real_en, end_site)

            # calcualte R1 start site  based on SE mode and strand
            start_site = utr.Strand == "+" ? r1_pos["start"] : r1_pos["end"]

            push!(st_arr, start_site - utr_site)
            push!(real_st, start_site)
        end

        if length(st_arr) > min_reads
           return TestData(utr = utr, st_arr = st_arr, en_arr = en_arr, real_st = real_st, real_en = real_en)
        end

        return nothing
    end

    function get_record_from_bam(
        path::Union{Vector, String}, 
        utr;
        single_end_mode::Bool = false,
        min_reads::Int=0,
        cage_mode::Bool = false
    )::Union{TestData, ExtractData, Nothing}

        if cage_mode
            return read_from_bam_cage(path, utr, min_reads=min_reads)
        end
        return read_from_bam(path, utr, single_end_mode=single_end_mode, min_reads=min_reads)
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


    function is_reads_pass(exon_coord::Dict, exon_range::AbstractArray, site_relative::AbstractArray; is_read1_junc::Bool, is_read2_junc::Bool)::Bool
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

        return pass1 && pass2 && !failed_on_border
    end

    function get_record_from_bam_transcript(path::String, transcript, exons::Vector; min_reads::Int=0, expand::String="200,1000")::Union{TestData, Nothing}

        sites = Dict{}

        try
            expand = parse.(Int, split(expand, ","))
            if length(expand) < 2
                throw(string("Only single expand value was passed: ", expand, "; 200,1000 format was required"))
            end
        catch e
            throw(e)
        end

        st_arr, en_arr = Vector(), Vector()
        real_st, real_en = Vector(), Vector()
        r1, r2 = Dict(), Dict()

        utr = Genomic.BED(
            transcript.Chrom,
            transcript.Strand == "+" ? max(1, transcript.Start - expand[1]) : max(1, transcript.Start - expand[2]),
            transcript.End == "+" ? transcript.End + expand[2] : transcript.End + expand[1],
            transcript.Name, transcript.GeneID,
            transcript.Strand
        )
    
        utr_site = utr.Strand == "+" ? utr.Start : utr.End

        bai = string(path, ".bai")
        reader = open(BAM.Reader, path, index=bai)
        for record in eachoverlap(reader, utr.Chrom, utr.Start:utr.End)
            if !filter(record)
                continue
            end

            # Only kept R1
            if BAM.flag(record) & SAM.FLAG_READ2 > 0
                r2[BAM.tempname(record)] = Dict(
                    "start" => BAM.leftposition(record), 
                    "end" => BAM.rightposition(record),
                    "junc" => occursin("N", BAM.cigar(record))
                )
                continue
            end

            if utr.Start > BAM.leftposition(record) || utr.End < BAM.rightposition(record)
                continue
            end

            # read strand
            strand = determine_strand(record)

            if strand != utr.Strand
                continue
            end

            start_site = utr.Strand == "+" ? BAM.leftposition(record) : BAM.rightposition(record)

            r1[BAM.tempname(record)] = Dict(
                "start"=>BAM.leftposition(record), 
                "end"=>BAM.rightposition(record),
                "junc" => occursin("N", BAM.cigar(record))
            )
        end
        close(reader)


        # convert genomic site of exons to relative pos in transcript
        exon_coord = Dict{Int, Int}()
        exon_range = Vector()
        for e = exons
            for i = e.Start:e.End
                exon_coord[i] = exons[1].Start + length(exon_coord) + 1
            end
            append!(exon_range, [e.Start, e.End])
        end
        
        utr_site = utr.Strand == "+" ? utr.Start : utr.End
        for (name, r1_pos) = r1
            r2_pos = get(r2, name, nothing)

            if isnothing(r2_pos)
                continue
            end

            # R1 needs locate in UTR
            if utr.Start > r1_pos["start"] || utr.End < r1_pos["end"]
                continue
            end

            site_relative = [
                convert_absolute_relative(exon_coord, exon_range, [r1_pos["start"], r1_pos["end"]], r1_pos["junc"], is_read2 = false)...
                convert_absolute_relative(exon_coord, exon_range, [r2_pos["start"], r2_pos["end"]], r2_pos["junc"], is_read2 = true)...
            ]

            if !any(isnan.(site_relative))
                  
                if !(r1_pos["junc"] || r2_pos["junc"]) || is_reads_pass(exon_coord, exon_range, site_relative, is_read1_junc=r1_pos["junc"], is_read2_junc=r2_pos["junc"])
                    # R2 end site, for Rn
                    end_site = strand == "-" ? site_relative[3] : site_relative[4]
                    # calcualte R1 start site  based on SE mode and strand
                    start_site = utr.Strand == "+" ? site_relative[1] : site_relative[2]
                    
                    push!(st_arr, start_site - utr_site)
                    push!(real_st, strand == "+" ? r1_pos["start"] : r1_pos["end"])

                    push!(en_arr, end_site - utr_site)
                    push!(real_en, strand == "-" ? r2_pos["start"] : r2_pos["end"])
                    continue
                end
            end
        end

        if length(st_arr) > min_reads
            return TestData(utr = utr, st_arr = st_arr, en_arr = en_arr, real_st = real_st, real_en = real_en, exon_coord = exon_coord)
        end
        
        return nothing
    end
end
