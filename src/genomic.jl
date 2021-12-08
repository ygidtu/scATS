import Base.isless


module Genomic
    using CodecZlib
    using Formatting
    using Parameters
    using ProgressMeter

    export GTF, BED, load_UTRs, new

    @with_kw struct GTF
        Chrom::String
        Start::Int64
        End::Int64
        Strand::String
        Type::String
        ID::String
        Name::String
        GeneID::String
        TranscriptID::String
        Attributes::Dict{String, String}
    end

    Base.show(io::IO, g::GTF) = print(io,  g.Chrom, ":", g.Start, "-", g.End, ":", g.Strand, "\t", g.Type, "\t", g.ID, "-", g.Name, "\t", g.GeneID, "\t", g.TranscriptID)

    # GTF related functions
    # convert gtf attributes to dict
    function extract_attributes(record::SubString)::Dict{String, String}
        res = Dict()
        temp = split(record, " ")

        for i = 1:2:length(temp)
            if i+1 <= length(temp)
                res[temp[i]] = strip(temp[i+1], ['"', ';'])
            end
        end

        for k = ["gene_name", "transcript_name"]
            if !haskey(res, k)
                res[k] = get(res, "Parent", "")
            end
        end

        for k = ["gene_id", "transcript_id"]
            if !haskey(res, k)
                res[k] = get(res, "Parent", "")
            end
        end

        return res
    end
    
    # create gtf from gtf line
    function create_GTF(record::String)::GTF

        records = split(strip(record), "\t")
        attributes = extract_attributes(records[9])
        
        id, name, gene_id, transcript_id = "", "", "", ""
        if records[3] == "gene"
            id, name = attributes[string(records[3], "_id")], get(attributes, string(records[3], "_name"), "")
        elseif records[3] in ["transcript", "exon"]
            id, name, gene_id, transcript_id = attributes[string(records[3], "_id")], get(attributes, string(records[3], "_name"), ""), attributes["gene_id"], attributes["transcript_id"]
        end

        return GTF(
            records[1], 
            parse(Int64, records[4]),
            parse(Int64, records[5]), 
            records[7],
            records[3], id, name, gene_id, transcript_id,
            attributes
        )
    end

    function load_UTRs(gtf::AbstractString; utr_length::Int = 500)::Vector
        data = Dict(
            "transcript" => Dict{String, Genomic.GTF}(),  # transcript id and GTF
            "exon" => Dict{String, Vector}()  # transcript id and it's chiledren exons
        )

        open(gtf) do io
            seekend(io)
            fileSize = position(io)
            seekstart(io)

            p = Progress(fileSize, 1, "Loading...")   # minimum update interval: 1 second

            stream = eachline(io)
            if endswith(gtf, ".gz")
                stream = eachline(GzipDecompressorStream(io))
            end

            for line = eachline(io)
                if startswith(line, "#")
                    continue
                end

                gtf = Genomic.create_GTF(line)

                if occursin(r"(transcript|RNA)"i, gtf.Type)
                    data[gtf.Type][gtf.ID] = gtf
                elseif occursin(r"(exon)"i, gtf.Type)
                    if !haskey(data["transcript"], gtf.TranscriptID)
                        continue
                    end

                    if !haskey(data[gtf.Type], gtf.TranscriptID)
                        data[gtf.Type][gtf.TranscriptID] = Vector()
                    end
                    push!(data[gtf.Type][gtf.TranscriptID], gtf)
                end
                update!(p, position(io))
            end
        end

        # calculate UTR by exons range
        utr_length = div(utr_length, 2)
        utrs = Dict()

        l = Threads.SpinLock()
        p = Progress(length(data["transcript"]), 1, "Calculating UTRs...")
        Threads.@threads for transcript = collect(values(data["transcript"]))
            lock(l)
            if !haskey(utrs, transcript.Strand)
                utrs[transcript.Strand] = Vector{BED}()
            end
            unlock(l)

            exons = sort(get(data["exon"], transcript.ID, []))

            if length(exons) < 1
                continue
            end

            site = exons[1].Start
            if transcript.Strand != "+"
                site = exons[length(exons)].End 
            end

            lock(l)
            push!(utrs[transcript.Strand], BED(
                Chrom = transcript.Chrom,
                Start = site - utr_length,
                End = site + utr_length,
                Name = transcript.Attributes["gene_name"],
                Score = transcript.Name,
                Strand = transcript.Strand
            ))
            unlock(l)

            next!(p)
        end

        # merging
        res = Vector()
        Threads.@threads for strand = collect(keys(utrs))
            temp = merge(utrs[strand])

            lock(l)
            res = vcat(res, temp)
            unlock(l)
        end

        return sort(res)
    end

    @with_kw struct BED
        Chrom::String
        Start::Int64
        End::Int64
        Name::String
        Score::String
        Strand::String
        Sites::Vector{Int} = []
    end

    Base.show(io::IO, b::BED) = print(io, join([b.Chrom, b.Start, b.End, b.Name, b.Score, b.Strand], "\t"))
    Base.:(==)(x::BED, y::BED) = x.Chrom == y.Chrom && x.Start == y.Start && x.End == y.End && x.Strand == y.Strand

    function new(
        chrom::AbstractString, start_pos::Int, end_pos::Int, strand::AbstractString; 
        name::AbstractString=".", score::AbstractString=".", sites::Vector = []
    )
        return BED(string(chrom), start_pos, end_pos, string(name), string(score), string(strand), sites)
    end

    # merge BED
    function merge(beds::Vector{BED})::Vector{BED}
        res = Vector()

        beds = sort(beds)

        if length(beds) < 2
            return beds
        end

        old_bed = beds[1]
        for i = 2:length(beds)
            if isupstream(old_bed, beds[i])
                push!(res, old_bed)
                old_bed = beds[i]
            elseif isdownstream(old_bed, beds[i])
                push!(res, beds[i])
            else
                old_bed = BED(
                    Chrom = old_bed.Chrom,
                    Start = old_bed.Start, 
                    End = beds[i].End,
                    Name = join(unique([old_bed.Name, beds[i].Name]), ","), 
                    Score = join(unique([old_bed.Score, beds[i].Score]), ","),
                    Strand = old_bed.Strand
                )
            end
        end
        push!(res, old_bed)
        return res
    end

    function unique_beds(beds::Vector)::Vector
        res = Vector()

        beds = sort(beds)

        if length(beds) <= 1
            return beds
        end

        curr = beds[1]
        for i in 2:length(beds)
            if beds[i] != curr
                push!(res, curr)
                curr = beds[i]
            end
        end

        push!(res, curr)

        return res
    end

    Region = Union{Genomic.GTF, Genomic.BED}

    # Region related functions
    function isupstream(a::Region, b::Region)::Bool
        if a.Chrom != b.Chrom
            return a.Chrom < b.Chrom
        end

        return a.End < b.Start
    end

    function isdownstream(a::Region, b::Region)::Bool
        if a.Chrom != b.Chrom
            return a.Chrom > b.Chrom
        end

        return a.Start > b.End
    end

    function isoverlap(a::Region, b::Region)::Bool
        return a.Chrom == b.Chrom && a.Start < b.End && a.End > b.Start
    end
end


function isless(a::Genomic.Region, b::Genomic.Region)
    if a.Chrom != b.Chrom
        return a.Chrom < b.Chrom
    end

    if a.Start != b.Start
        return a.Start < b.Start
    end

    if a.End != b.End
        return a.End < b.End
    end
    return a.Strand < b.Strand
end

