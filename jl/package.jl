using Pkg

packages = ["Memento", "KernelEstimator", "StatsBase", "ArgParse", "FilePathsBase", "ProgressMeter", "Formatting", "ProgressBars", "Distributions", "RCall"]

for i in sort(packages)
    println(i)
    Pkg.add(i)
end

precompile
# julia -e 'using Pkg; Pkg.add("Memento"); Pkg.add("")'