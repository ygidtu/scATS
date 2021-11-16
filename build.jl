#!/usr/bin/env julia

#=
This script is used to install julia requirements
=#

using Pkg


function main()
    prepare = false
    deps = Vector{String}()
    open(string(@__DIR__, "/Project.toml")) do r
        while !eof(r)
            line = readline(r)
            if startswith(line, "[deps]")
                prepare = true
                continue
            end

            if prepare
                push!(deps, strip(split(line, "=")[1]))
            end
        end
        close(r)
    end

    println(string("install deps: ", deps))

    Pkg.add(deps)

    println("precompile")
    precompile
end

main()