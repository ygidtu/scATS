module Bam
    using BioAlignments
    using Compat: @__MODULE__  # requires a minimum of Compat 0.26. Not required on Julia 0.7
    using FilePathsBase
    using Formatting
    using GenomicFeatures
    using Libz
    using Memento
    using Parameters
    using ProgressMeter
    using XAM


    include(joinpath(@__DIR__, "genomic.jl"))

    # Create our module level LOGGER (this will get precompiled)
    const LOGGER = Memento.config!("info"; fmt="[{date} - {level} | {name}]: {msg}")
    # const LOGGER = getlogger(@__MODULE__)

    # Register the module level LOGGER at runtime so that folks can access the LOGGER via `get_LOGGER(MyModule)`
    # NOTE: If this line is not included then the precompiled `MyModule.LOGGER` won't be registered at runtime.
    function __init__()
        Memento.register(LOGGER)
    end

    function setLevel(level::String)
        setlevel!(LOGGER, level)
    end

    @with_kw struct BamData
        path::String
        alias::Union{String, Nothing} = nothing
        barcodes::Dict{String, Set{String}} = Dict()
    end

    Base.show(io::IO, self::BamData) = print(
        io,
        Formatting.format(
            FormatExpr("{}\t{}"), # \nreal_st: {}\nreal_en: {}
            self.alias, self.path
        )
    )

    function prepare_bam_list(path::String)::Vector{BamData}
        fs = Vector{BamData}()
        for line in eachline(open(path))
            line = split(line, "\t")

            if !exists(Path(line[1]))
                error(LOGGER, Formatting.format(
                    FormatExpr("{} not exists"), line[1]
                ))
            end

            alias = ""
            if length(line) > 1
                alias = line[2]
            end

            bamPath = absolute(Path(line[1]))
            if islink(line[1])
                bamPath = readlink(bamPath)
            end

            push!(fs, BamData(string(bamPath), alias, load_barcodes(bamPath)))
        end

        return fs
    end

    function load_barcodes(path::AbstractPath)::Dict{String, Set}
        barcodes = Dict()
        path = joinpath(parent(path), "filtered_feature_bc_matrix/barcodes.tsv.gz")

        if !exists(path)
            return barcodes
        end

        for line in eachline(open(path) |>  ZlibInflateInputStream)
            line = strip(line)
            key = line[1:min(length(line), 3)]

            if !haskey(barcodes, key)
                barcodes[key] = Set()
            end

            push!(barcodes[key], line)
        end
        return barcodes
    end

    function check_barcode(bam::BamData, barcode::AbstractString)::Bool
        if length(bam.barcodes) < 1
            return true
        end

        key = barcode[1:min(3, length(barcode))]

        return haskey(bam.barcodes, key) && (barcode in bam.barcodes[key])
    end

    function filter(record, bulk::Bool = false)::Bool
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
    
        if BAM.flag(record) & SAM.FLAG_PROPER_PAIR == 0
            return false
        end
        
        if bulk
            return true
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
    
    function count_reads(bam_list::Vector, region, cell_tag::String="CB", umi_tag::String="UB")::Dict{String, Int}
        res = Dict() 
    
        for bam in bam_list
            reader = open(BAM.Reader, bam.path, index=string(bam.path, ".bai"))
            for record in eachoverlap(reader, region.Chrom, region.Start:region.End)
                if !filter(record)
                    continue
                end
        
                # read strand
                strand = determine_strand(record)
        
                if strand != region.Strand
                    continue
                end
        
                if BAM.leftposition(record) < region.End && BAM.rightposition(record) > region.Start
                    auxdata = Dict(BAM.auxdata(record))
                    cb = auxdata[cell_tag]
                    ub = auxdata[umi_tag]
    
                    if check_barcode(bam, cb)
                        if bam.alias != ""
                            cb = string(bam.alias, "_", cb)
                        end

                        if !haskey(res, cb)
                            res[cb] = Dict()
                        end
            
                        res[cb][ub] = get(res[cb], ub, 0) + 1
                    end
                end
            end
        end
        return Dict(x => length(y) for (x, y) = res)
    end

    function count_bulk(bam_list::Vector, region)::Dict{String, Int}

        res = Dict()
        for bam in bam_list
            reader = open(BAM.Reader, bam.path, index=string(bam.path, ".bai"))
            count = 0
            for record in eachoverlap(reader, region.Chrom, region.Start:region.End)
                if !filter(record)
                    continue
                end
        
                # read strand
                strand = determine_strand(record)
        
                if strand != region.Strand
                    continue
                end
        
                if BAM.leftposition(record) < region.End && BAM.rightposition(record) > region.Start
                    count += 1
                end
            end
            res[bam.alias] = count
        end
        return res
    end

    # function test()
    #     for i in prepare_bam_list("/mnt/raid61/Personal_data/zhangyiming/code/afe/tests/covid.txt")
    #         println(i)
        
    #         for (key, value) in reads(i, Genomic.BED("1", 182460, 182461, "test", "1", "+"))
    #             println(key)
    #             println(value)
    #             break
    #         end
    #         break
    #     end
    # end

end

# Bam.test()

