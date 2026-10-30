# Shared trace checks include independent small-system exact enumerators.
include("../multimodel/verify.jl")
include("common.jl")
function verify_scaling(root;require_current_source=false)
    verify(root;require_current_source)
    config=TOML.parsefile(joinpath(root,"config.toml"))
    manifest=TOML.parsefile(joinpath(root,"manifest.toml"))
    @testset "critical design and connectivity provenance" begin
        @test haskey(manifest,"finished_utc")
        @test file_hash(joinpath(root,"config.toml"))==manifest["config_sha256"]
        cases=critical_cases(config)
        @test length(cases)==length(manifest["runs"])
        for (index,(case,row)) in enumerate(zip(cases,manifest["runs"]))
            @test all(row[string(k)]==v for (k,v) in pairs(case))
            @test row["seed"]==config["seed"]+1009index
            @test row["p"]==sqrt(case.q)/(1+sqrt(case.q))
            @test row["beta"]==log1p(sqrt(case.q))
            meta=TOML.parsefile(joinpath(root,"runs",row["run_id"],"metadata.toml"))
            @test meta["schema_version"]==3
            @test meta["connectivity_backend"]==(case.family=="RandomCluster" ? "interleaved_bfs" : "not_applicable")
        end
        @test length(manifest["restart_audits"])==3
        for audit in manifest["restart_audits"]
            @test file_hash(joinpath(root,"checkpoints",audit["run_id"]*".checkpoint"))==audit["checkpoint_sha256"]
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    paths=filter(!=("--require-current-source"),ARGS)
    root=isempty(paths) ? joinpath(@__DIR__,"results","critical") : abspath(first(paths))
    verify_scaling(root;require_current_source="--require-current-source" in ARGS)
end
