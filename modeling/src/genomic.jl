import Base.isless


module Genomic
    using BSON
    using Parameters
    using ProgressMeter

    export GTF, Bed, create_GTF, load_GTF, isupstream, isdownstream, get_bed, new_bed

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

    @with_kw struct BED
        Chrom::String
        Start::Int64
        End::Int64
        Name::String
        Score::String
        Strand::String
    end

    Base.show(io::IO, b::BED) = print(io, join([b.Chrom, b.Start, b.End, b.Name, b.Score, b.Strand], "\t"))

    function get_bed(bed::BED, expand::Int=100)::String
        return join([
            bed.Chrom, 
            max(0, bed.Start - expand), 
            bed.End + expand, 
            bed.Name, bed.Score, 
            bed.Strand
        ], "\t")
    end

    function new_bed(line::String)::BED  # ; chromosomes::Dict{String, String}=nothing
        lines = split(strip(line), "\t")
        
        if length(lines) >= 6
            return BED(
                replace(lines[1], "chr"=>""), parse(Int64, lines[2]), parse(Int64, lines[3]),
                lines[4], strip(lines[5]), lines[6]
            )
        elseif length(lines) == 3
            return BED(
                lines[1], parse(Int64, lines[2]), parse(Int64, lines[3]),
                ".", ".", "."
            )
        elseif length(lines) == 4
            return BED(
                lines[1], parse(Int64, lines[2]), parse(Int64, lines[3]),
                ".", ".", strip(lines[4])
            )
        else
            error(LOGGER, string("the number of columns mismatch: ", lines))
            exit(1)
        end
    end

    function get_bed_short(self::BED)::String
        return Formatting.format(
            FormatExpr("{}:{}-{}:{}"),
            self.chrom, self.start_pos, self.end_pos, self.strand
        )
    end

    Region = Union{Genomic.GTF, Genomic.BED}

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

    # sort gtf
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
                    Name = join(unique([old_bed.Name, beds[i].Name]), "|"), 
                    Score = join(unique([old_bed.Score, beds[i].Score]), "|"),
                    Strand = old_bed.Strand
                )
            end
        end
        push!(res, old_bed)
        return res
    end

    function load_GTF(gtf::AbstractString)::Dict
        data = Dict(
            "transcript"=>Dict{String, Genomic.GTF}(), 
            "exon"=>Dict{String, Vector}()
        )

        open(gtf) do io
            seekend(io)
            fileSize = position(io)
            seekstart(io)

            p = Progress(fileSize, 1)   # minimum update interval: 1 second
            while !eof(io)
                line = readline(io)
                if startswith(line, "#")
                    continue
                end

                gtf = Genomic.create_GTF(line)

                if gtf.Type == "transcript"
                    # if there is not gene_biotype or gene_biotype is not specific types passed
                    if !haskey(gtf.Attributes, "gene_biotype") || !(gtf.Attributes["gene_biotype"] in ["antisense", "lincRNA", "protein_coding"])
                        continue
                    end

                    data[gtf.Type][gtf.ID] = gtf
                elseif gtf.Type == "exon"
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

        return data
    end

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

    return a.End < b.End
end

