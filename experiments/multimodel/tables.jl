"""Readers for this study's explicit, unquoted numeric CSV schema.

This is deliberately not a general CSV implementation. Reject duplicate column
names, ragged rows and unsupported cells instead of silently coercing data.
"""
module MultiModelTables
export read_table, write_table, read_trace

function read_table(path)
    open(path) do io
        names = String.(split(readline(io), ','))
        length(unique(names)) == length(names) || error("duplicate CSV columns: $path")
        records = Dict{String,String}[]
        for line in eachline(io)
            cells = split(line, ',')
            length(cells) == length(names) || error("ragged CSV row: $path")
            push!(records, Dict(name => String(cell) for (name, cell) in zip(names, cells)))
        end
        records
    end
end

function write_table(path, records, fields)
    open(path, "w") do io
        println(io, join(fields, ','))
        for record in records
            cells = [string(record[field]) for field in fields]
            any(cell -> occursin(',', cell) || occursin('\n', cell), cells) &&
                error("unsupported unquoted CSV cell")
            println(io, join(cells, ','))
        end
    end
end

function read_trace(path)
    rows = read_table(path)
    length(rows) > 1 || error("trace needs at least two observations")
    result = Dict(name => [parse(Float64, row[name]) for row in rows]
                  for name in keys(first(rows)))
    result["measurement"] == collect(1:length(rows)) || error("invalid measurement ordering")
    all(values -> all(isfinite, values), values(result)) || error("nonfinite trace")
    result
end
end
