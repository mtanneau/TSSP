import Base.*

mutable struct Realization
    p::Float64
    values::Dict{Tuple{String,String},Float64}
end

*(r1::Realization, r2::Realization) = Realization(r1.p*r2.p, merge(r1.values, r2.values))

function parse_sto_file(file_name)
    name = ""
    block_name = ""
    block_realization = []

    indeps = Dict{Tuple{String, String}, Vector{Realization}}()
    blocks = Dict{String, Vector{Realization}}()

    realization = Realization(0.0, Dict{Tuple{String,String},Float64}())

    open(file_name) do f

        section = :SEC

        for ln in eachline(f)
            s = split(ln, r" +")

            if s[1] == "STOCH"
                section = :STOCH
                name = s[2]
                continue
            elseif s[1] == "INDEP"
                # INDEP section
                section = :INDEP
                continue
            elseif s[1] == "BLOCKS"
                section = :BLOCKS
                continue
            elseif s[1] == "ENDATA"
                # End of file reached
                break
            else
                s[1] == "" || error("Unexpected line:\n$(ln)\n$(s)")
            end

            if section == :INDEP
                parse_indep_line!(s, indeps)
            end

            if section == :BLOCKS
                if s[2] == "BL"
                    # New block entry
                    # Get block name and probability
                    block_name = s[3]
                    
                    block_prob = parse(Float64, s[5])
                    realization = Realization(block_prob, Dict{Tuple{String,String},Float64}())

                    if !haskey(blocks, block_name)
                        blocks[block_name] = Realization[]
                    end

                    push!(blocks[block_name], realization)
                else
                    # Read values
                    parse_block_line!(s, realization)
                end
            end

        end
    end

    return name, indeps, blocks
end

# Parse INDEP section
function parse_indep_line!(s, indeps)

    if length(s) == 5
        # Special case for `fxm2_6.sto` and `fxm2_16.sto`
        colname = String(s[2])
        rowname = String(s[3])
        coeffval = parse(Float64, s[4])
        prob = parse(Float64, s[5])

    elseif length(s) >= 6
        # Normal case
        colname = String(s[2])
        rowname = String(s[3])
        coeffval = parse(Float64, s[4])
        prob = parse(Float64, s[6])
    end

    r = Realization(prob, Dict((rowname, colname) => coeffval))

    if !haskey(indeps, (colname, rowname))
        indeps[colname, rowname] = Realization[]
    end

    push!(indeps[colname, rowname], r)

    return indeps
end


function parse_block_line!(s, block_realization)

    colname = String(s[2])
    rowname = String(s[3])
    val = parse(Float64, s[4])
    block_realization.values[rowname, colname] = val

    # optional: extra row
    if length(s) == 6
        rowname = String(s[5])
        val = parse(Float64, s[6])
        block_realization.values[rowname, colname] = val
    end

    return block_realization
end