module StartPSI
    using FilePathsBase
    using Formatting
    using Base
    export new, run, Data

    struct Data
        Input::String
        Gtf::String
        Junction::String
        Output::String
        Length::Int64
        ExonicGff::String
        ExonicInclusion::String
        ExonicExclusion::String
        FilteredJunction::String
        PSI::String
    end

    function new(input::String, gtf::String, junction::String, output::String, len::Int64)::Data

        return Data(
            input, gtf, junction, output, len,
            string(output, "_exonic_parts.gff"),
            string(output, "_exonic_parts.inclusion"),
            string(output, "_exonic_parts.exclusion"),
            string(output, "_filtered_junctions.bed"),
            string(output, "_exonic_parts.psi")
        )

    end

    function CountInclusion(self::Data)
        # fe = "{}\t{}:{}\n"
        # open(self.ExonicGff, "w+") do w
        #     open(self.Gtf, "r") do io
        #         while !eof(io)
        #             lines = split(replace(strip(readline(io)), r"[\";]" => ""), "\t")
        #             write(w, Formatting.format(fe, join(lines[1:8], "\t"), lines[14], lines[12]))
        #         end
        #         close(io)
        #     end
        #     close(w)
        # end

        reads = self.Input
        gff = self.Gtf

        open(pipeline(`coverageBed -nonamecheck -split -abam $reads -b $gff`, `awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$5-$4+1,$9,$10}'`, `sort -k5,5`), "r", stdout) do io
            open(self.ExonicInclusion, "w+")  do w
                while !eof(io)
                    write(w, readline(io))
                end
                close(w)
            end
        end
    end

    function FilterJunctions(self::Data)
        open(self.Junction, "r") do io
            open(self.FilteredJunction, "w+") do w
                while !eof(io)
                    line = strip(readline(io))
                    if !occursin("description", line)
                        lines = split(replace(line, "," => "\t"), "\t")
                        write(w, join([
                            lines[1],
                            string(parse(Int64, lines[2])+parse(Int64, lines[13])), 
                            string(parse(Int64, lines[3])-parse(Int64, lines[14])),
                            lines[4], lines[5], lines[6]
                        ], "\t"))
                        write(w, "\n")
                    end
                end
                close(w)
            end
            close(io)
        end
    end

    function CountExclusion(self::Data)
        gff = self.Gtf
        junc = self.FilteredJunction
        open(self.ExonicExclusion, "w+") do w
            open(pipeline(`bedtools intersect -nonamecheck -wao -f 1.0 -s -a $gff -b $junc`, `awk 'BEGIN{OFS="\t"}{$16 == 0? s[$9] += 0:s[$9] += $14}END{for (i in s) {print i,s[i]}}'`, `sort -k1,1`), "r", stdout) do io
                while !eof(io)
                    write(w, readline(io))
                end
            end
            close(w)
        end
        rm(self.FilteredJunction)
    end

    function start(self::Data)
        FilterJunctions(self)
        CountInclusion(self)
        CountExclusion(self)

        # curr_dir = @__DIR__
        # bash = string(curr_dir, "/../ref/ExonicPartPSI_2.sh")

        # GFF=self.Gtf
        # INBAM=self.Input
        # readLength=self.Length
        # JUNCTIONS=self.Junction
        # PREFIX=self.Output
        # Base.run(`bash $bash bedtools $GFF $INBAM $readLength $JUNCTIONS $PREFIX`)
    end
end