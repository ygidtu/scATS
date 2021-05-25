#!/usr/bin/env julia
#=
Created at 2021.05.21
=#


module Isoform
    include(joinpath(@__DIR__, "struct.jl"))

    # function map_iso_to_gene(iso_list::WinList.Data, gene_list::WinList.Data)

    #     if !iso_list.isSorted
    #         WinList.sort(iso_list)
    #     end

    #     if !gene_list.isSorted
    #         WinList.sort(gene_list)
    #     end

    #     res_win_list = WinList.new()
    #     curr_exon_j = 1
    #     local_iso_index = 1

    #     for win = gene_list.Data
    #         if curr_exon_j >= length(iso_list.Data)
    #             res_win_list = WinList.append(res_win_list, Window.new())
    #             continue
    #         end

    #         if Window.isLessEqual(win, iso_list.Data[curr_exon_j])
    #             res_win_list = WinList.append(res_win_list, Window.shiftStart(win, local_iso_index))
    #             local_iso_index += len(win.Data)

    #             if iso_list.Data[curr_exon_j].End == win.End
    #                 curr_exon_j += 1
    #             end
    #         else
    #             res_win_list = append(res_win_list, Window.new())
    #         end
    #     end

    #     return res_win_list
    # end

end

include(joinpath(@__DIR__, "struct.jl"))

using .DataStruct




function main()

    iso_list = WinListUtils.create()

    for i in 100:200:1000
        iso_list = WinListUtils.append(iso_list, WindowUtils.create(i, i + 100))
    end

    gene_list = WinListUtils.create()
    for i in 100:300:1000
        iso_list = WinListUtils.append(gene_list, WindowUtils.create(i, i + 200))
    end

    println(iso_list)
    println(gene_list)

    iso_list = WinListUtils.rmDuplicate(iso_list)
    println(iso_list)

end


main()
