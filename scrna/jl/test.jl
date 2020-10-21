include("src/utils.jl")


function main()
    f1 = "d_exon.bed"
    a = read(`bedtools subtract -a $f1 -b d_intron.bed`, String)
    a = split(a, "\n")
    println(string("test ", a[1]))
end

# main()
println(gen_k_arr(10, 2))