

module DataStruct
    using Parameters

    @with_kw mutable struct Window
        Start::Int64
        End::Int64
    end

    Base.show(io::IO, win::Window) = print(io,  win.Start, "\t", win.End, "\t", win.End - win.Start)

    function WindowStr(win::Window)::String
        return string(win.Start, "|", win.End, "|", win.End - win.Start)
    end

    function isless(a::Window, b::Window)::Bool
        if a.Start == b.Start
            return a.End < b.End
        end

        return a.Start < b.Start
    end


    @with_kw mutable struct WinList
        Data::Vector{Window}
        isSorted::Bool
    end


    @with_kw struct TreeNode
        Left
        Right
        WinList::WinList
     end

    export Window, WinList, WindowUtils, WinListUtils, TreeNode


    module WindowUtils
        using ..DataStruct:Window

        export Window

        function isEmpty(self::Window)::Bool
            return self.Start == self.End
        end

        function Len(self::Window)::Int64
            return self.End - self.Start
        end

        function lShift(self::Window, other::Window)::Bool
            return self.End <= other.Start
        end        

        function isLess(self::Window, other::Window)::Bool
            return self.Start < other.Start < self.End < other.End
        end

        function isLessEqual(self::Window, other::Window)::Bool
            return other.Start <= self.Start && self.End <= other.End
        end

        function isEqual(self::Window, other::Window)::Bool
            return self.Start == other.Start && self.End == other.End
        end

        function isGreatEqual(self::Window, other::Window)::Bool
            return self.Start <= other.Start && other.End <= self.End
        end

        function isGreat(self::Window, other::Window)::Bool
            return other.Start < self.Start < other.End < self.End
        end

        function rShift(self::Window, other::Window)::Bool
            return other.End <= self.Start
        end

        function notEqual(self::Window, other::Window)::Bool
            return ! (other.End <= self.Start || self.End <= other.Start)
        end

        function adjacent(self::Window, other::Window)::Bool
            return self.End == other.Start || self.Start == other.End
        end

        function contains(self::Window, point::Int64)::Bool
            if isEmpty(self)
                return false
            end

            return self.Start <= point < self.End
        end

        function shiftStart(self::Window, start::Int64)::Window
            if self.Start == 1 << 31
                return self
            end

            return Window(start, start + self.End - self.Start)
        end

        function create(startSite::Int64, endSite::Int64)::Window
            return Window(startSite, endSite)
        end

        function new()::Window
            return Window(1 << 31, 1 << 31)
        end
    end


    module WinListUtils

        using ..DataStruct:Window, WinList, WindowStr, WindowUtils, isless

        function create()
            return WinList([], false)
        end

        function append(data::WinList, win::Window)
            push!(data.Data, win)
            data.isSorted = false
            return data
        end

        function rmEmpty(data::WinList)
            data.Data = [x for x = data.Data if !WindowUtils.isEmpty(x::Window)]
            return data
        end

        function sort(data::WinList)
            data.Data = Base.sort(data.Data, lt=isless)
            return data
        end

        function rmDuplicate(data::WinList)
            temp = Dict{String, Window}()

            for i = data.Data
                temp[WindowStr(i)] = i
            end

            data.Data = [x for x = values(temp)]
            return sort(data)
        end

        function merge(self::WinList)
            winlist = rmEmpty(self)

            if ! winlist.isSorted
                winlist = sort(winlist)
            end

            resList = append(create(), winlist[1])

            for win = winlist
                curr_win = resList[length(resList)]

                if curr_win.Start <= win.Start <= curr_win.End
                    resList[length(resList)] = WindowUtils.create(curr_win.Start, win.End)
                else
                    resList = append(resList, win)
                end
            end

            resList.isSorted = true

            return resList
        end

        function split(self::WinList)
            winlist = rmEmpty(self)

            if length(winlist.WinList) == 0
                return winlist
            end

            if !winlist.isSorted
                winlist = sort(winlist)
            end

            boarder = Dict{Int64, String}

            for win = winlist
                boarder[win.Start] = ""
                boarder[win.End] = ""
            end

            boarder_arr = Base.sort(keys(border))

            winlist = create()
            for i = 1:(length(boarder_arr) - 1)
                winlist = append(winlist, WindowUtils.create(boarder_arr[i], boarder_arr[i+1]))
            end
            winlist.isSorted = true
            return winlist
        end

        function getLeft(self::WinList)
            for win = self.WinList
                if win.Start != 1 << 31
                    return win.Start
                end
            end
            return 1 << 31
        end


        function getRight(self::WinList)
            for win = reversed(self.WinList)
                if win.End != 1 << 31
                    return win.End
                end
            end
            return 1 << 31
        end

        function getRange(self::WinList)
            return WindowUtils.create(getLeft(self), getRight(self))
        end

    end

end