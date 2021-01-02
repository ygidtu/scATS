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
    using XAM

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
    
    function new_bed(line::String; chromosomes::Dict{String, String}=nothing)::BEDRecord
        lines = split(strip(line), ",")
        chrom = lines[2]
        if !isnothing(chromosomes) && length(chromosomes) > 0
            if !haskey(chromosomes, chrom)
                if startswith(chrom, "chr")
                    chrom = replace(chrom, "chr"=>"")
                else
                    chrom = string("chr", chrom)
                end

                if !haskey(chromosomes, chrom)
                    return BEDRecord("", 0, 0, "", "", "")
                end
            end
        end

        return BEDRecord(
            chrom, parse(Int64, lines[3]), parse(Int64, lines[4]),
            lines[6], lines[11], lines[1]
        )
        
        # if length(lines) >= 6
        #     return BEDRecord(
        #         lines[1], parse(Int64, lines[2]), parse(Int64, lines[3]),
        #         lines[6], lines[4], lines[5]
        #     )
        # elseif length(lines) == 3
        #     return BEDRecord(
        #         lines[1], parse(Int64, lines[2]), parse(Int64, lines[3]),
        #         "", "", ""
        #     )
        # elseif length(lines) == 4
        #     return BEDRecord(
        #         lines[1], parse(Int64, lines[2]), parse(Int64, lines[3]),
        #         strip(lines[4]), "", ""
        #     )
        # else
        #     error(LOGGER, string("the number of columns mismatch: ", lines))
        #     exit(1)
        # end
    end

    struct ExtractedData
        utr::BEDRecord
        site1::String
        site2::String
        r1_utr_st_arr::Vector
        r1_len_arr::Vector
    end

    function get_hash(data::ExtractedData)::String
        return string(join(sort(data.r1_utr_st_arr), ","), "|", join(sort(data.r1_len_arr), ","))
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
        path::String, promoters::Array;
        max_number_of_sites::Int64=-1, seed::Int64=42,
        distance::Int64=-500, min_reads::Int=10,
        # min_distance::Int=10, max_distance::Int=500
    )::Vector{ExtractedData}
        # reads_start, reads_end = Vector{String}(), Vector{String}()
        # start_sites = Dict{Int64, Int}()
    
        bai = string(path, ".bai")
        # println(string(utr, " start"))

        promoter_sites = Dict{BEDRecord, Vector}()
        reader = open(BAM.Reader, path, index=bai)

        sites = Dict{}
        for utr = promoters
            
        end


        res = Vector()
        for i = 1:length(promoters)
            for j = (i+1):length(promoters)

                # if start_sites_keys[j] - start_sites_keys[i] < min_distance || start_sites_keys[j] - start_sites_keys[i] > min_distance
                #     continue
                # end

                # if start_sites[start_sites_keys[i]] < min_reads && start_sites[start_sites_keys[i]] < min_reads
                #     continue
                # end
                
                # r1_utr_st_arr = [x[1] for x = data if x[3] == start_sites_keys[i] || x[3] == start_sites_keys[j]]
                # r1_len_arr = [x[2] for x = data if x[3] == start_sites_keys[i] || x[3] == start_sites_keys[j]]
                
                start_site = promoters[i].start_pos
                if  promoters[i].strand == "-"
                    start_site = promoters[i].end_pos
                end

                end_site = promoters[j].start_pos
                if  promoters[j].strand == "-"
                    start_site = promoters[j].end_pos
                end

                utr = BEDRecord(
                    promoters[i].chrom,
                    start_site, end_site,
                    promoters[i].strand,
                    string(promoters[i].score, ",", promoters[j].score),
                    promoters[i].name
                )

                r1_utr_st_arr = []
                r1_len_arr = []
                
                for record in eachoverlap(reader, utr.chrom, utr.start_pos:utr.start_pos)
                    if !filter(record)
                        continue
                    end
                    
                    # Only kept R2
                    if BAM.flag(record) & SAM.FLAG_READ1 > 0
                        continue
                    end

                    if BAM.leftposition(record) < start_site && BAM.rightposition(record) > end_site
                        continue
                    end

                    utr_site = utr.start_pos
                    start_site = BAM.position(record)
                    if utr.strand != "+"
                        start_site = BAM.rightposition(record)
                    end

                    if BAM.position(record) <= utr_site <= BAM.rightposition(record) && start_site - utr_site > distance
                        # pos in UTR
                        push!(r1_utr_st_arr, start_site - utr_site)

                        # length in UTR
                        push!(r1_len_arr, min(BAM.rightposition(record), utr.end_pos) - max(BAM.position(record), utr.start_pos)) 
                    end
                end

                if length(r1_utr_st_arr) <= min_reads
                    continue
                end

                if max_number_of_sites > 0 && length(r1_len_arr) > max_number_of_sites
                    r1_len_arr = randsubseq(Random.seed!(seed), r1_len_arr, max_number_of_sites / length(r1_len_arr))
                    r1_utr_st_arr = randsubseq(Random.seed!(seed), r1_utr_st_arr, max_number_of_sites / length(r1_utr_st_arr))
                end
                
                push!(
                    res, 
                    ExtractedData(
                        utr, 
                        promoters[i].name,
                         promoters[j].name, 
                        r1_utr_st_arr, r1_len_arr,
                    )
                )
            end
        end

        # for promoter = promoters
        #     data = Vector{Vector{Int64}}()
        #     try
        #         for record in eachoverlap(reader, promoter.chrom, promoter.start_pos:promoter.end_pos)
        #             if !filter(record)
        #                 continue
        #             end
                    
        #             # Only kept R2
        #             if BAM.flag(record) & SAM.FLAG_READ1 > 0
        #                 continue
        #             end

        #             strand = promoter.strand

        #             utr_site = promoter.end_pos
        #             start_site = BAM.position(record)
        #             if strand != "+"
        #                 utr_site = promoter.start_pos
        #                 start_site = BAM.rightposition(record)
        #             end

        #             # start_sites[start_site] = get(start_sites, start_site, 0) + 1

        #             if BAM.position(record) <= utr_site <= BAM.rightposition(record) && start_site - utr_site > distance
        #                 t1 = start_site-utr_site                                                                      # pos in UTR
        #                 t2 = min(BAM.rightposition(record), promoter.end_pos) - max(BAM.position(record), promoter.start_pos)   # length in UTR

        #                 push!(data, [t1, t2, start_site])
        #                 # push!(reads_start, string(BAM.position(record)))
        #                 # push!(reads_end, string(BAM.rightposition(record)))
        #             end
        #         end
                
        #         promoter_sites[promoter] = data
        #     catch e
        #         warn(LOGGER, string(promoter))
        #         error(LOGGER, e)
        #         return ""
        #     end
        # end
        close(reader)

        # data = unique(data)
        # maxmium used 200 (default) sites to save memory
        # if length(data) > max_number_of_sites
        #     Random.seed!(seed)
        #     data = Random.shuffle(data)
        #     data = data[1:max_number_of_sites]
        # end
        
        #=
        if length(promoter_sites)  < 2
            return Vector()
        end

        
        if length(promoter_sites) > 1
            # start_sites_keys = sort(collect(keys(start_sites)))
        
            for i = 1:length(promoters)
                for j = (i+1):length(promoters)

                    # if start_sites_keys[j] - start_sites_keys[i] < min_distance || start_sites_keys[j] - start_sites_keys[i] > min_distance
                    #     continue
                    # end

                    # if start_sites[start_sites_keys[i]] < min_reads && start_sites[start_sites_keys[i]] < min_reads
                    #     continue
                    # end
                    
                    # r1_utr_st_arr = [x[1] for x = data if x[3] == start_sites_keys[i] || x[3] == start_sites_keys[j]]
                    # r1_len_arr = [x[2] for x = data if x[3] == start_sites_keys[i] || x[3] == start_sites_keys[j]]

                    r1_utr_st_arr = []
                    r1_len_arr = []
                    for x = promoter_sites[promoters[i]]
                        push!(r1_utr_st_arr, x[1])
                        push!(r1_len_arr, x[2])
                    end
                    for x = promoter_sites[promoters[j]]
                        push!(r1_utr_st_arr, x[1])
                        push!(r1_len_arr, x[2])
                    end
                    
                    if length(r1_utr_st_arr) < 2
                        continue
                    end

                    push!(
                        res, 
                        ExtractedData(
                            # utr, 
                            promoters[i].name, promoters[j].name, 
                            abs(promoters[j].end_pos - promoters[i].start_pos),
                            r1_utr_st_arr, r1_len_arr,
                            # reads_start, reads_end
                        )
                    )
                end
            end
        end
        =#

        return res
    end

    function run(
        data::Vector{ExtractedData},
        mu_f::Int64=300, min_ws::Float64=0.01, 
        min_pa_gap::Int64=100, max_beta::Int64=70, theta_step::Int64=9, 
        sigma_f::Int64=50, using_R::Bool=false, verbose::Bool=true
    )::Vector{String}
        res = Vector{String}()
        if length(data) < 0
            return res
        end

        try
            # utr = data[1].utr
            r1_len_arr = data[1].r1_len_arr
            r1_utr_st_arr  = data[1].r1_utr_st_arr

            self = APA.new(
                5, 1,
                r1_utr_st_arr .+ abs(data[1].utr.end_pos - data[1].utr.start_pos),
                r1_len_arr,
                [NaN for i = 1:length(r1_utr_st_arr)],
                [NaN for i = 1:length(r1_utr_st_arr)],
                [NaN for i = 1:length(r1_utr_st_arr)],
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

            if isnothing(temp_res.bic)
                return res
            end

            for d = data
                push!(
                    res, 
                    Formatting.format(
                        FormatExpr("{}:{}-{}:{}\t{}-{}\t{}"), 
                        d.utr.chrom, d.utr.start_pos, d.utr.end_pos, d.utr.strand,
                        d.site1, d.site2,
                        # join(d.reads_starts, ","), join(d.reads_ends, ","), 
                        string(temp_res)
                    )
                )
            end
        catch e
            info(LOGGER, string(e))
            if verbose
                debug(LOGGER, string(e))
            end
        end
        return res
    end
end