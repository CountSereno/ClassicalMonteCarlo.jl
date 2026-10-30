using ClassicalMonteCarlo, Random, Dates
include("common.jl")

function specification(case,config,seed;burn=config["burn"],samples=config["samples"])
    common=(;L=case.L,initial=Symbol(case.initial),algorithm=Symbol(case.algorithm),seed,burn,samples,stride=config["stride"])
    case.family=="Potts" ? PottsRunSpec(;common...,q=Int(case.q),beta=log1p(sqrt(case.q))) :
        RandomClusterRunSpec(;common...,q=case.q,p=sqrt(case.q)/(1+sqrt(case.q)))
end
physical_state(model::BondConfiguration)=model.occupied
physical_state(model::Potts)=model.spins
function audit_restart(original,output,id)
    partial=start_run(original.spec)
    midpoint=original.spec.burn+original.spec.samples÷2*original.spec.stride
    advance!(partial,midpoint)
    path=joinpath(output,"checkpoints",id*".checkpoint")
    save_checkpoint(path,partial)
    resumed=run!(load_checkpoint(path))
    equal=observations(resumed)==observations(original) &&
        physical_state(resumed.model)==physical_state(original.model) &&
        resumed.measurement_updates==original.measurement_updates &&
        resumed.kernel_work==original.kernel_work && resumed.pending_work==original.pending_work &&
        resumed.updates==original.updates && rand(copy(resumed.rng),UInt64,16)==rand(copy(original.rng),UInt64,16)
    equal || error("restart trajectory mismatch")
    Dict("run_id"=>id,"exactly_equal"=>equal,"checkpoint_update"=>midpoint,
         "checkpoint_sha256"=>file_hash(path),"scope"=>"observations, final physical configuration, schedule/work counters and future RNG")
end

function run_study(config_path,output)
    config=TOML.parsefile(config_path)
    ispath(output) && error("choose a new output directory")
    all(b->config["samples"]%b==0 && config["samples"]÷b>=2,config["block_sizes"]) || error("invalid blocks")
    hashes=generation_sources()
    cases=critical_cases(config)
    for family in ("RandomCluster","Potts")
        case=first(filter(c->c.family==family,cases))
        run!(start_run(specification(case,config,0;burn=4,samples=8)))
    end
    mkpath(output); cp(config_path,joinpath(output,"config.toml"))
    manifest=Dict{String,Any}("schema_version"=>1,"started_utc"=>string(now(UTC)),
        "julia_version"=>string(VERSION),"machine"=>Sys.MACHINE,"cpu_name"=>Sys.CPU_NAME,
        "threads"=>Threads.nthreads(),"boundary"=>"periodic square","source_sha256"=>hashes,
        "config_sha256"=>file_hash(config_path),"runs"=>Any[],"restart_audits"=>Any[],
        "timing"=>"sequential production plus measurements; compilation, burn, explicit pre-timer GC and I/O excluded")
    for (index,case) in enumerate(cases)
        id="run_"*lpad(index,3,'0'); seed=config["seed"]+1009index
        state=start_run(specification(case,config,seed))
        advance!(state,state.spec.burn)
        GC.gc()
        seconds=@elapsed run!(state)
        export_run(joinpath(output,"runs",id),state)
        row=Dict{String,Any}(string(k)=>v for (k,v) in pairs(case))
        merge!(row,Dict("run_id"=>id,"seed"=>seed,"samples"=>config["samples"],"burn"=>config["burn"],
            "stride"=>config["stride"],"p"=>sqrt(case.q)/(1+sqrt(case.q)),"beta"=>log1p(sqrt(case.q)),
            "nvertices"=>case.L^2,"nedges"=>2case.L^2,"production_seconds"=>seconds,
            "update_unit"=>update_unit(state.kernel),"work_unit"=>work_unit(state.kernel)))
        push!(manifest["runs"],row)
        if index==1 || (case.L==maximum(config["sizes"]) && case.q==3 && case.chain==1)
            push!(manifest["restart_audits"],audit_restart(state,output,id))
        end
        manifest["measurements"]=sum(r["samples"] for r in manifest["runs"])
        save_toml(joinpath(output,"manifest.toml"),manifest)
        println("[$index/$(length(cases))] $id $(case.family) q=$(case.q) L=$(case.L): $(round(seconds;digits=3)) s")
        flush(stdout)
    end
    hashes==generation_sources() || error("source changed during production")
    manifest["finished_utc"]=string(now(UTC))
    save_toml(joinpath(output,"manifest.toml"),manifest)
end

if abspath(PROGRAM_FILE)==@__FILE__
    config=isempty(ARGS) ? joinpath(@__DIR__,"configs","critical.toml") : abspath(ARGS[1])
    output=length(ARGS)<2 ? joinpath(@__DIR__,"results","critical") : abspath(ARGS[2])
    run_study(config,output)
end
