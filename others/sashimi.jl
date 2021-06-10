using ArgParse
using RCall
using DataFrames
using FilePathsBase


struct BED
    Chrom::String
    Start::Int
    End::Int
    Strand::String
end


function create(data::AbstractString)::BED
    sites = split(data, ":")

    pos = split(sites[2], "-")

    chrom = split(sites[1], "-")
    chrom = chrom[length(chrom)]

    if occursin("|", chrom)
        chrom = split(chrom, "|")
        chrom = chrom[length(chrom)]
    end

    return BED(
        chrom,
        parse(Int, pos[1]),
        parse(Int, pos[2]),
        sites[3]
    )
end


function toRegion(data::BED; expand::Int=1000)::String
    return string(data.Chrom, ":", data.Start - expand, ":", data.End + expand)
end


struct Sites
    sites::Vector
    gene::String
    region::String
end


function create_site(line::String, expand::Int = 1000)
    try
        utr, site, gene, n = split(line, ",")
        utr = create(utr)
        site = create(site)
    
        return Sites(
            [site.Strand == "+" ? site.Start : site.End],
            gene, toRegion(utr, expand=expand)
        )
    catch e
        utr, site, gene = split(line, "\t")
        utr = create(utr)
        site = create(site)
    
        return Sites(
            [site.Strand == "+" ? site.Start : site.End],
            gene, toRegion(utr, expand=expand)
        )
    end

end


function concat_sites(s1::Sites, s2::Sites)::Sites
    sites = [s1.sites..., s2.sites...]

    return Sites(sites, s1.gene, s1.region)
end


function load_filtered_list(path::String)::Dict
    res = Dict()
    idx = Dict()
    open(path) do r
        while !eof(r)
            line = readline(r)

            site = create_site(line)

            line = split(line, ",")
            if length(line) < 2
                line = split(line[1], "\t")
            end
            idx[line[2]] = line[1]

            if !haskey(res, line[1])
                res[line[1]] = site
            else
                res[line[1]] = concat_sites(res[line[1]], site)
            end
        end
    end
    return Dict(x => res[y] for (x, y) = idx)
end


function make_sashimi(data::Sites; gtf::String, bam::String, bc::String, order::String, output::String, expand::Int=1000)

    # region = toRegion(data, expand=expand)

    # line = data.Strand == "+" ? data.Start : data.End

    pdf = joinpath(output, string(data.gene, "_", replace(data.region, ":" => "_"), ".pdf"))
    region = data.region
    line = join(map(string, data.sites), ",")

    if exists(Path(pdf))
        return
    end
    
    try
        run(`sashimiplot junc --gtf $gtf --bam $bam --bc $bc --sj 100 --co $order --junc $region --fileout $pdf --trackline $line --ie 1,1  --ps RF --ssm R1`)
    catch e
        println(e)
        println(`sashimiplot junc --gtf $gtf --bam $bam --bc $bc --sj 100 --co $order --junc $region --fileout $pdf --trackline $line --ie 1,1  --ps RF --ssm R1`)
    end
    
end


#=
dist = round(Int, site.End - site.Start)

# lines = join([string(site.Start), string(site.End)], ",")
lines = site.Site
event = string(site.Chrom, ":", site.Start - dist, ":", site.End + dist) # , ":", site.Strand
focus = string(site.Start, "-", site.End)
dots = site.RealST
dots2 = site.RealEn
run(`sashimiplot junc --gtf $ref --bam $bam --sj 1000 --junc $event --ie 1,1  --ps RF --ssm R1 --fileout $o --trackline $lines --focus $focus`)


gene = obj@assays$RNA@counts["HLA-DPB1", rownames(obj@meta.data)[obj@meta.data$CellType == "CD14 mono"]]

gene = as.data.frame(gene)
gene$stage = obj@meta.data[rownames(gene), "recovery_stage"]
gene$log2 = log2(gene$gene + 1)
gene$log10 = log10(gene$gene + 1)

p1 = ggplot(gene, aes(x=stage, y=gene, fill=stage)) +
    geom_violin()

p2 = ggplot(gene, aes(x=stage, y=log2, fill=stage)) +
    geom_violin()  

p3 = ggplot(gene, aes(x=stage, y=log10, fill=stage)) +
    geom_violin()

p4 = ggplot(gene, aes(x=stage, y=gene, fill=stage)) +
    geom_violin() +
    scale_y_log10()

ggsave(
    filename = "/mnt/raid64/Covid19_Gravida/HLA-DPB1_vln.pdf",
    plot = cowplot::plot_grid(p1, p2, p3, p4, ncol = 2),
    width = 8, height = 8
)
=#


function main(rds::String, gtf::String, bam::String, bc::String, order::String, output::String; expand::Int=1000, cellType::Union{Nothing, String}=nothing, n::Int=10)
    markers = rcopy(R"""
        dat = readRDS($rds)
        if (!'gene' %in% colnames(dat)) {
            dat$gene = rownames(dat)
            dat = dat[, c('gene', 'Tau')]
        }
        dat
    """)

    if "Tau" in names(markers)
        markers = markers[markers.Tau .> .8, :]
    elseif "CellType" in names(markers)
        markers = markers[markers.p_val_adj .< 0.05, :]
        markers = markers[[!(x in ["platelet", "erythroid", "HSC"]) for x = markers.CellType], :]
    end

    temp = markers

    if "CellType" in names(markers)
        if !isnothing(cellType) && cellType != ""
            markers = markers[markers.CellType .== replace(cellType, "_"=>" "), :]
        end

        temp = combine(groupby(markers, :CellType)) do sdf
            col = "avg_logFC" in names(markers) ? :avg_logFC : :deltaPSI
            first(sort(sdf, col, rev="deltaPSI" in names(markers)), n)
        end
    end

    # if "comp" in names(markers)
    #     markers = markers[occursin.("HC", markers.comp), :]
    #     temp = combine(groupby(markers, :comp)) do sdf
    #         col = "avg_logFC" in names(markers) ? :avg_logFC : :deltaPSI
    #         first(sort(sdf, col, rev="deltaPSI" in names(markers)), n)
    #     end
    # end
    try
        if !exists(Path(output))
            mkdir(Path(output), recursive=true)
        end
    catch e
        error(logger, Formatting.format(FormatExpr("Error while create {}: {}"), output, e))
        exit(1)
    end

    sites = load_filtered_list("/mnt/raid64/ATS/Personal/zhangyiming/overall_ATS_nslc.postprocess")  # "/mnt/raid61/Personal_data/zhangyiming/code/afe/others/filtered_list.csv"
    println(length(sites))
    println(size(temp))
    Threads.@threads for row in eachrow(temp) # Threads.@threads 
        key = split(row.gene, "|")
        key = key[2]
        try
            if haskey(sites, key)
                make_sashimi(sites[key], gtf=gtf, bam=bam, bc=bc, order=order, output=output)
            end
        catch e
            println(row)
            println(e)
        end
    end

end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "Path to markers rds"
            arg_type = String
            required = true
        "--gtf", "-g"
            help = "Path to gtf file"
            arg_type = String
            required = true
        "--bam", "-b"
            help = "Path to bam file"
            arg_type = String
            required = true
        "--bc", "-c"
            help = "Path to barcode list file"
            arg_type = String
            required = true
        "--order", "-r"
            help = "Path to barcode list file"
            arg_type = String
            required = true
        "--cell"
            help = "Path to barcode list file"
            arg_type = String
        "--expand", "-e"
            help = "region to expand around the ATS"
            arg_type = Int
            default = 1000
        "--top", "-n"
            help = "region to expand around the ATS"
            arg_type = Int
            default = 10
        "--output", "-o"
            help = "Prefix of output file"
            arg_type = String
            required = true
    end

    return parse_args(s)
end

args = parse_commandline()

main(
    get(args, "input", ""), 
    get(args, "gtf", ""),
    get(args, "bam", ""),
    get(args, "bc", ""),
    get(args, "order", ""),
    get(args, "output", ""),
    expand=get(args, "expand", 1000),
    cellType = get(args, "cell", nothing),
    n = get(args, "top", 10)
)


#=



sort(data, :Tau, rev = true)

=#