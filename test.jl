using BenchmarkTools


function TestPmap()
    
    res = pmap(1:10000) do i
        return i
    end

    open("test.txt", "w+") do w
        for r = res
            write(w, string(r, "\n"))
        end
        close(w)
    end
end


function TestPmap1()

    open("test1.txt", "w+") do w
        Threads.@threads for i = 1:10000
            write(w, string(i, "\n"))
        end
        close(w)
    end
end

@btime TestPmap1()