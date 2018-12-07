function parse_time_file(file_name)
    name = ""
    is_lp = false
    n_periods = 1

    col_periods = String[]  # Name of beginning column for each period
    row_periods = String[]  # Name of beginning row for each period

    open(file_name) do f
        for ln in eachline(f)
            s = split(ln, r" +")
            if s[1] == "TIME"
                name = s[2]
            end

            if s[1] == ""
                # read name
                push!(col_periods, s[2])
                push!(row_periods, s[3])
            end

            if s[1] == "ENDATA"
                # End of file reached
                break
            end
        end
    end

    return name, is_lp, col_periods, row_periods
end