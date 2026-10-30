# Read-only artifact checks, plus a closed-form weighted-fit reference that does
# not call production fitting code. Run after analysis and plotting.
using Test, Statistics
include("common.jl")
include("../multimodel/tables.jl")
using .MultiModelTables

function verify_analysis(root;require_current_source=false)
    report=TOML.parsefile(joinpath(root,"analysis","report.toml"))
    figures=TOML.parsefile(joinpath(root,"figures","metadata.toml"))
    number(row,key)=parse(Float64,row[key])
    @testset "saved inputs, tables, figures and source fingerprints" begin
        for (path,hash) in report["input_sha256"]
            @test file_hash(joinpath(root,path))==hash
        end
        for (path,hash) in report["output_sha256"]
            @test file_hash(joinpath(root,"analysis",path))==hash
        end
        for (path,hash) in figures["input_sha256"]
            @test file_hash(joinpath(root,path))==hash
        end
        for (path,hash) in figures["output_sha256"]
            @test file_hash(joinpath(root,"figures",path))==hash
        end
        @test length(figures["output_sha256"])==12
        if require_current_source
            for document in (report,figures), (path,hash) in document["source_sha256"]
                @test file_hash(joinpath(REPOSITORY,path))==hash
            end
            @test file_hash(joinpath(REPOSITORY,"benchmarks/results/connectivity/summary.csv"))==figures["benchmark_sha256"]
        end
    end
    aggregate=read_table(joinpath(root,"analysis","aggregate.csv"))
    fits=read_table(joinpath(root,"analysis","power_fits.csv"))
    @testset "independent weighted slope and error reference" begin
        for fit in fits
            selected=filter(row->row["family"]=="RandomCluster" && row["q"]==fit["q"] &&
                row["observable"]==fit["observable"] && number(row,"L")>=number(fit,"Lmin"),aggregate)
            x=[log(number(row,"L")) for row in selected]
            y=[log(number(row,"estimate")) for row in selected]
            w=[(number(row,"estimate")/number(row,"stderr"))^2 for row in selected]
            xbar,ybar=sum(w.*x)/sum(w),sum(w.*y)/sum(w)
            denominator=sum(w.*(x.-xbar).^2)
            exponent=sum(w.*(x.-xbar).*(y.-ybar))/denominator
            @test exponent ≈ number(fit,"exponent") atol=1e-12
            @test inv(sqrt(denominator)) ≈ number(fit,"stderr") atol=1e-12
            @test number(fit,"dof")==length(selected)-2
            @test number(fit,"bootstrap_low")<=number(fit,"bootstrap_high")
        end
    end
    println("Verified derived artifacts and independent weighted fits.")
end
if abspath(PROGRAM_FILE)==@__FILE__
    paths=filter(!=("--require-current-source"),ARGS)
    root=isempty(paths) ? joinpath(@__DIR__,"results","critical") : abspath(first(paths))
    verify_analysis(root;require_current_source="--require-current-source" in ARGS)
end
