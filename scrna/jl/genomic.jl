import Base.isless


module Genomic
    export GTF, Bed, create_GTF, sort_gtf

    struct GTF
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

    Base.show(io::IO, g::GTF) = print(
        io, 
        g.Chrom, ":", g.Start, "-", g.End, ":", g.Strand,
        "\t", g.Type, "\t", g.ID, "-", g.Name, "\t", g.GeneID, "\t", g.TranscriptID
    )

    struct BED
        Chrom::String
        Start::Int64
        End::Int64
        Name::String
        Score::String
        Strand::String
    end

    Base.show(io::IO, b::BED) = print(
        io, 
        join([b.Chrom, b.Start, b.End, b.Name, b.Score, b.Strand], "\t")
    )

    # convert gtf attributes to dict
    function extract_attributes(record::SubString)::Dict{String, String}
        res = Dict()
        temp = split(record, " ")

        for i = 1:2:length(temp)
            if i+1 <= length(temp)
                res[temp[i]] = strip(temp[i+1], ['"', ';'])
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
            id, name = attributes[string(records[3], "_id")], attributes[string(records[3], "_name")]
        elseif records[3] in ["transcript", "exon"]
            id, name, gene_id, transcript_id = attributes[string(records[3], "_id")], get(attributes, string(records[3], "_name"), ""), attributes["gene_id"], attributes["transcript_id"]
        end

        return GTF(
            records[1], parse(Int64, records[4]),
            parse(Int64, records[5]), records[7],
            records[3], id, name, gene_id, transcript_id,
            attributes
        )
    end

    # sort gtf
    function merge(beds::Vector{BED}, distance::Int64=0)::Vector{BED}
        res = Vector()

        beds = sort(beds)

        if length(beds) < 2
            return beds
        end

        old_bed = beds[1]
        for i = 2:length(beds)
            if old_bed < beds[i]
                push!(res, old_bed)
                old_bed = beds[i]
            elseif old_bed.End + distance > beds[i].Start
                old_bed = BED(
                    old_bed.Chrom,
                    old_bed.Start, beds[i].End,
                    old_bed.Name, old_bed.Score,
                    old_bed.Strand
                )
            else
                push!(res, old_bed)
                old_bed = beds[i]
            end
        end
        push!(res, old_bed)
        return res
    end
end

Region = Union{Genomic.GTF, Genomic.BED}

function isless(a::Region, b::Region)
    if a.Chrom != b.Chrom
        return a.Chrom < b.Chrom
    end

    if a.Start != b.Start
        return a.Start < b.Start
    end

    return a.End < b.End
end
