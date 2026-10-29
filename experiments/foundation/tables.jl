"""Small readers for the deliberately simple numeric CSV schema of this study.

These are not general CSV parsers: the generator emits unquoted cells without
commas/newlines. Reject ragged rows instead of silently misreading a trace.
"""
module FoundationTables
export read_table, write_table, read_trace

function read_table(path)
    open(path) do io
        names = split(readline(io), ',')
        records = Dict{String,String}[]
        for line in eachline(io)
            cells = split(line, ',')
            length(cells) == length(names) || error("ragged CSV row in $path")
            push!(records,Dict(String(k)=>String(v) for (k,v) in zip(names,cells)))
        end
        records
    end
end

function write_table(path,records,fields)
    open(path,"w") do io
        println(io,join(fields,','))
        for record in records
            cells = string.([record[f] for f in fields])
            any(s->occursin(',',s)||occursin('\n',s),cells) && error("unsupported CSV cell")
            println(io,join(cells,','))
        end
    end
end

function read_trace(path)
    records = read_table(path)
    n = length(records)
    n > 1 || error("trace requires at least two measurements")
    measurement = [parse(Int,r["measurement"]) for r in records]
    measurement == collect(1:n) || error("trace samples are incomplete or unordered")
    (E=[parse(Float64,r["total_energy"]) for r in records],
     M=[parse(Int,r["total_magnetization"]) for r in records],
     updates=[parse(Int,r["update"]) for r in records],
     work=[parse(Int,r["kernel_work"]) for r in records])
end
end
