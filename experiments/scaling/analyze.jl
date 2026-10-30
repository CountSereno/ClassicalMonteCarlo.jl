using ClassicalMonteCarlo, Statistics, Random, Printf
include("common.jl")
include("../multimodel/tables.jl")
using .MultiModelTables

record(;kwargs...)=Dict{String,Any}(string(k)=>v for (k,v) in kwargs)
theory(q,observable)=observable=="C1" ? (q==2 ? 15/8 : 28/15) : (q==2 ? 7/4 : 26/15)

function study_observables(row,trace)
    n=row["nvertices"]
    if row["family"]=="RandomCluster"
        Dict("C1"=>trace["largest_cluster"],
             "S2_over_N"=>trace["cluster_size_squared_sum"]./n,
             "energy_density"=>-trace["occupied_bonds"]./(row["p"]*n))
    else
        Dict("S2_over_N"=>n.*trace["simplex_squared"],
             "energy_density"=>trace["total_energy"]./n)
    end
end

function correlation_summary(values,maxlag)
    rho=autocorrelation(values;maxlag)
    window=0
    for lag in 1:maxlag
        rho[lag+1]>0 || break
        window=lag
    end
    tau=.5+sum(@view rho[2:window+1])
    (;rho,tau,window)
end

# Stratified block bootstrap: resample within each independent chain, preserving
# pairing within a block. Here each fit uses one observable across independent
# sizes. Quantiles describe statistical uncertainty conditional on the fit model.
function bootstrap_slopes(rng,sizes,means,errors,chain_blocks,replicas)
    slopes=Float64[]
    for _ in 1:replicas
        resampled=Float64[]
        for blocks in chain_blocks
            push!(resampled,mean(mean(rand(rng,b,length(b))) for b in blocks))
        end
        # Hold the original log-space weights fixed during resampling.
        fit=power_law_fit(sizes,resampled;sigma=resampled.*(errors./means))
        push!(slopes,fit.exponent)
    end
    quantile(slopes,[.025,.975])
end

function analyze(root)
    manifest=TOML.parsefile(joinpath(root,"manifest.toml"))
    haskey(manifest,"finished_utc") || error("production is incomplete")
    config=TOML.parsefile(joinpath(root,"config.toml"))
    destination=joinpath(root,"analysis"); mkpath(destination)
    primary=config["primary_block"]
    chain_rows=Dict{String,Any}[]; block_rows=Dict{String,Any}[]; acf_rows=Dict{String,Any}[]
    blocks=Dict{Tuple{String,Float64,Int,String},Vector{Vector{Float64}}}()
    discarded_blocks=Dict{Tuple{String,Float64,Int,String,Int},Vector{Vector{Float64}}}()
    input_hashes=Dict("manifest.toml"=>file_hash(joinpath(root,"manifest.toml")),"config.toml"=>file_hash(joinpath(root,"config.toml")))
    for row in manifest["runs"]
        relative=joinpath("runs",row["run_id"],"trace.csv")
        meta_relative=joinpath("runs",row["run_id"],"metadata.toml")
        input_hashes[relative]=file_hash(joinpath(root,relative))
        input_hashes[meta_relative]=file_hash(joinpath(root,meta_relative))
        input_hashes[relative]==TOML.parsefile(joinpath(root,meta_relative))["trace_sha256"] || error("trace hash mismatch")
        trace=read_trace(joinpath(root,relative))
        for (observable,values) in sort!(collect(study_observables(row,trace));by=first)
            correlation=correlation_summary(values,config["maxlag"])
            primary_blocks=block_means(values,primary)
            key=(row["family"],row["q"],row["L"],observable)
            push!(get!(blocks,key,Vector{Float64}[]),primary_blocks)
            for discard in (0,length(values)÷4,length(values)÷2)
                retained=block_means(values[discard+1:end],primary)
                push!(get!(discarded_blocks,(key...,discard),Vector{Float64}[]),retained)
            end
            half=length(primary_blocks)÷2
            early,late=primary_blocks[1:half],primary_blocks[half+1:end]
            drift_se=sqrt(var(early)/length(early)+var(late)/length(late))
            drift_z=(mean(early)-mean(late))/drift_se
            errors=Float64[]
            for blocksize in config["block_sizes"]
                b=block_means(values,blocksize)
                error=std(b)/sqrt(length(b)); push!(errors,error)
                push!(block_rows,record(;run_id=row["run_id"],family=row["family"],q=row["q"],L=row["L"],
                    observable,blocksize,stderr=error,blocks=length(b)))
            end
            stderr=std(primary_blocks)/sqrt(length(primary_blocks))
            push!(chain_rows,record(;run_id=row["run_id"],family=row["family"],q=row["q"],L=row["L"],
                chain=row["chain"],initial=row["initial"],observable,estimate=mean(values),stderr,
                tau=correlation.tau,window=correlation.window,ess=length(values)/(2correlation.tau),
                block_to_tau=primary/correlation.tau,large_block_ratio=last(errors)/stderr,
                drift_z,production_seconds=row["production_seconds"]))
            if row["L"]==maximum(config["sizes"])
                for lag in 0:config["maxlag"]
                    push!(acf_rows,record(;run_id=row["run_id"],family=row["family"],q=row["q"],L=row["L"],
                        observable,lag,acf=correlation.rho[lag+1]))
                end
            end
        end
        println("analyzed $(row["run_id"])"); flush(stdout)
    end
    aggregate=Dict{String,Any}[]
    for key in sort!(collect(keys(blocks)))
        family,q,L,observable=key
        chain_blocks=blocks[key]; chain_means=mean.(chain_blocks)
        n=length(chain_blocks)
        within=sqrt(sum(var(b)/length(b) for b in chain_blocks))/n
        between=std(chain_means)/sqrt(n)
        relevant=filter(r->(r["family"],r["q"],r["L"],r["observable"])==key,chain_rows)
        push!(aggregate,record(;family,q,L,observable,estimate=mean(chain_means),
            stderr=max(within,between),within_stderr=within,between_stderr=between,chains=n,
            mean_tau=mean(r["tau"] for r in relevant),max_tau=maximum(r["tau"] for r in relevant)))
    end
    discarded_aggregate=Dict{String,Any}[]
    for key in sort!(collect(keys(discarded_blocks)))
        family,q,L,observable,discard=key
        chains=discarded_blocks[key]; means=mean.(chains); n=length(chains)
        within=sqrt(sum(var(b)/length(b) for b in chains))/n
        between=std(means)/sqrt(n)
        push!(discarded_aggregate,record(;family,q,L,observable,discard,
            estimate=mean(means),stderr=max(within,between)))
    end
    discarded_fits=Dict{String,Any}[]
    for q in config["q_values"], observable in ("C1","S2_over_N"), discard in (0,config["samples"]÷4,config["samples"]÷2)
        rows=sort!(filter(r->r["family"]=="RandomCluster" && r["q"]==q && r["observable"]==observable &&
            r["discard"]==discard && r["L"]>=16,discarded_aggregate);by=r->r["L"])
        fit=power_law_fit([r["L"] for r in rows],[r["estimate"] for r in rows];sigma=[r["stderr"] for r in rows])
        push!(discarded_fits,record(;q,observable,discard,Lmin=16,exponent=fit.exponent,
            stderr=fit.exponent_stderr,chi2=fit.chi2,dof=fit.dof))
    end
    fits=Dict{String,Any}[]; effective=Dict{String,Any}[]; corrections=Dict{String,Any}[]
    covariance_rows=Dict{String,Any}[]
    rng=Xoshiro(config["analysis_seed"])
    for q in config["q_values"], observable in ("C1","S2_over_N")
        rows=sort!(filter(r->r["family"]=="RandomCluster" && r["q"]==q && r["observable"]==observable,aggregate);by=r->r["L"])
        sizes=Float64[r["L"] for r in rows]
        means=[r["estimate"] for r in rows]; errors=[r["stderr"] for r in rows]
        exponent=theory(q,observable)
        eff=effective_exponents(sizes,means;sigma=errors)
        for i in eachindex(eff.exponents)
            push!(effective,record(;q,observable,L_low=sizes[i],L_high=sizes[i+1],L_geometric=eff.sizes[i],
                exponent=eff.exponents[i],stderr=eff.stderr[i],theory=exponent))
            for j in eachindex(eff.exponents)
                push!(covariance_rows,record(;q,observable,i,j,covariance=eff.covariance[i,j]))
            end
        end
        for first_index in 1:length(rows)-2
            select=first_index:length(rows)
            ls,ys,es=sizes[select],means[select],errors[select]
            fit=power_law_fit(ls,ys;sigma=es)
            chain_blocks=[blocks[("RandomCluster",Float64(q),Int(L),observable)] for L in ls]
            quantiles=bootstrap_slopes(rng,ls,ys,es,chain_blocks,config["bootstrap_replicas"])
            push!(fits,record(;q,observable,Lmin=first(ls),Lmax=last(ls),n_sizes=length(ls),
                exponent=fit.exponent,stderr=fit.exponent_stderr,theory=exponent,
                amplitude=fit.amplitude,chi2=fit.chi2,dof=fit.dof,
                bootstrap_low=quantiles[1],bootstrap_high=quantiles[2]))
            # Fixed-power correction sensitivity. These are deliberately varied
            # nuisance ansatzes, not measurements of an irrelevant exponent.
            for omega in (.8,1.2,4/3,2.0)
                scaled=ys./ls.^exponent; scaled_errors=es./ls.^exponent
                corrected=finite_size_fit(ls,scaled;correction_exponent=omega,sigma=scaled_errors)
                push!(corrections,record(;q,observable,Lmin=first(ls),omega,
                    amplitude=corrected.intercept,amplitude_stderr=corrected.intercept_stderr,
                    correction=corrected.slope,chi2=corrected.chi2,dof=corrected.dof))
            end
        end
    end
    block_fits=Dict{String,Any}[]
    for q in config["q_values"], observable in ("C1","S2_over_N"), blocksize in config["block_sizes"]
        rows=sort!(filter(r->r["family"]=="RandomCluster" && r["q"]==q && r["observable"]==observable && r["L"]>=16,aggregate);by=r->r["L"])
        errors=Float64[]
        for row in rows
            selected=filter(r->r["family"]=="RandomCluster" && r["q"]==q && r["observable"]==observable &&
                r["L"]==row["L"] && r["blocksize"]==blocksize,block_rows)
            within=sqrt(sum(r["stderr"]^2 for r in selected))/length(selected)
            push!(errors,max(within,row["between_stderr"]))
        end
        fit=power_law_fit([r["L"] for r in rows],[r["estimate"] for r in rows];sigma=errors)
        push!(block_fits,record(;q,observable,blocksize,Lmin=16,exponent=fit.exponent,stderr=fit.exponent_stderr,
            chi2=fit.chi2,dof=fit.dof))
    end
    cross=Dict{String,Any}[]
    for L in config["spin_sizes"], observable in ("energy_density","S2_over_N")
        a=only(filter(r->r["family"]=="Potts" && r["L"]==L && r["observable"]==observable,aggregate))
        b=only(filter(r->r["family"]=="RandomCluster" && r["q"]==config["spin_q"] && r["L"]==L && r["observable"]==observable,aggregate))
        push!(cross,record(;L,q=config["spin_q"],observable,spin=a["estimate"],spin_stderr=a["stderr"],
            bond=b["estimate"],bond_stderr=b["stderr"],z=(a["estimate"]-b["estimate"])/hypot(a["stderr"],b["stderr"])))
    end
    tables=("chains"=>chain_rows,"block_errors"=>block_rows,"acf"=>acf_rows,"aggregate"=>aggregate,
            "power_fits"=>fits,"effective_exponents"=>effective,"effective_covariance"=>covariance_rows,
            "correction_fits"=>corrections,"fk_comparison"=>cross,
            "discard_aggregate"=>discarded_aggregate,"discard_power_fits"=>discarded_fits,"block_power_fits"=>block_fits)
    for (name,rows) in tables
        write_table(joinpath(destination,name*".csv"),rows,sort!(collect(keys(first(rows)))))
    end
    report=Dict{String,Any}("measurements"=>manifest["measurements"],"chains"=>length(manifest["runs"]),
        "primary_block"=>primary,"bootstrap_replicas"=>config["bootstrap_replicas"],
        "max_abs_fk_z"=>maximum(abs(r["z"]) for r in cross),
        "max_tau"=>maximum(r["tau"] for r in chain_rows),
        "block_to_tau_below_10"=>count(r->r["block_to_tau"]<10,chain_rows),
        "large_block_flagged_chains"=>length(unique(r["run_id"] for r in chain_rows if r["large_block_ratio"]>1.5)),
        "large_block_ratio_above_1p5"=>count(r->r["large_block_ratio"]>1.5,chain_rows),
        "acf_window_at_limit"=>count(r->r["window"]==config["maxlag"],chain_rows),
        "half_drift_abs_z_above_3"=>count(r->abs(r["drift_z"])>3,chain_rows),
        "input_sha256"=>input_hashes,
        "source_sha256"=>source_hashes(["experiments/scaling/analyze.jl","experiments/scaling/common.jl",
            "experiments/multimodel/tables.jl","src/analysis.jl","src/scaling.jl"]),
        "output_sha256"=>Dict(name*".csv"=>file_hash(joinpath(destination,name*".csv")) for (name,_) in tables))
    save_toml(joinpath(destination,"report.toml"),report)
    write_report(root,manifest,config,report,fits,cross,block_fits)
    println("Analysis saved; max FK discrepancy $(report["max_abs_fk_z"]) combined SE")
end

function write_report(root,manifest,config,report,fits,cross,block_fits)
    open(joinpath(root,"RESULTS.md"),"w") do io
        println(io,"# Critical cluster scaling: executed results\n")
        println(io,"$(report["chains"]) independent-seed chains and $(report["measurements"]) retained measurements. Periodic square lattices, q=2 and 3, L=$(join(config["sizes"],", ")). Each chain: $(config["burn"]) warmup sweeps and $(config["samples"]) production records. Sampling at the known critical p=sqrt(q)/(1+sqrt(q)); no critical point is estimated.\n")
        println(io,"## Leading-power fits and size sensitivity\n")
        println(io,"C1 is the largest FK component size; S2_over_N is the sum of squared sizes of **all** components divided by N. The reference powers are respectively 15/8 and 7/4 (q=2), 28/15 and 26/15 (q=3). See [definitions and sources](../../README.md).\n")
        println(io,"| q | Observable | Lmin | Fitted power ± statistical SE | Reference | chi2 / dof |\n|---:|---|---:|---:|---:|---:|")
        for r in fits
            r["Lmin"] in (8,16,32) || continue
            @printf(io,"| %.0f | %s | %.0f | %.6f ± %.6f | %.6f | %.2f / %d |\n",r["q"],r["observable"],r["Lmin"],r["exponent"],r["stderr"],r["theory"],r["chi2"],r["dof"])
        end
        println(io,"\nFits use log-space weighted least squares with delta-method relative errors. Size errors are the larger of the combined within-chain block SE and between-chain SE (four chains). They are estimated, so chi2 is diagnostic. Adjacent effective exponents share data: their full covariance is saved. All size cutoffs, $(config["bootstrap_replicas"]) stratified block-bootstrap quantiles and fixed-power correction sensitivities (omega=0.8,1.2,4/3,2) are saved in analysis/. No correction exponent is inferred.\n")
        selected=filter(r->r["Lmin"]==16,fits)
        largest_reference_offset=maximum(abs(r["exponent"]-r["theory"])/r["stderr"] for r in selected)
        @printf(io,"For the illustrative L>=16 range, the four fitted powers are within %.2f diagnostic standard errors of their reference values. Raising Lmin widens errors because fewer sizes remain; the q=3 second-moment slope also shifts visibly between Lmin=8 and 16. A near-straight log-log curve alone therefore does not establish asymptotic scaling.\n\n",largest_reference_offset)
        println(io,"The fit errors and bootstrap intervals omit finite-size model bias and possible residual nonstationarity. Compare cutoff dependence, scaled amplitudes and residuals before interpreting a fitted slope as an asymptotic exponent. This study is a reproducible scaling demonstration, not a precision determination.\n")
        println(io,"## Independent Potts / FK checks\n\n| L | Observable | Spin ± SE | Bond ± SE | Difference / combined SE |\n|---:|---|---:|---:|---:|")
        for r in cross
            @printf(io,"| %d | %s | %.6f ± %.6f | %.6f ± %.6f | %.3f |\n",r["L"],r["observable"],r["spin"],r["spin_stderr"],r["bond"],r["bond_stderr"],r["z"])
        end
        println(io,"\nPotts SW and Sweeny are independently seeded. Energy uses -<b>/(pN); the second moment uses N*<simplex_squared>=<S2>/N. Maximum-colour order is never identified with the largest FK cluster.\n")
        println(io,"## Statistical diagnostics and reproducibility\n")
        println(io,"Primary blocks: $(report["primary_block"]) sweeps; alternatives: $(join(config["block_sizes"],", ")). Largest estimated tau: $(round(report["max_tau"];digits=2)) sweeps. Flags: block/tau<10: $(report["block_to_tau_below_10"]); largest-block SE/primary SE>1.5: $(report["large_block_ratio_above_1p5"]); ACF window exhausted: $(report["acf_window_at_limit"]); first-half/second-half discrepancy >3 SE: $(report["half_drift_abs_z_above_3"]). The large-block flags occur in $(report["large_block_flagged_chains"]) distinct chains. These are multiple, correlated diagnostics, not tests that prove equilibrium. Raw trajectories and every per-chain diagnostic are retained. `discard_power_fits.csv` repeats the L>=16 fits after discarding 0, one-quarter or one-half of production; these overlapping fits are correlated and do not replace longer burn tests.\n")
        println(io,"Reweighting the L>=16 fits with each block length gives the following sensitivity check. Means are unchanged; weights and error estimates change. At block length 2048 only eight blocks remain per chain, so those error estimates are themselves noisy. No flag is silently removed.\n")
        println(io,"| q | Observable | Block length | Power ± statistical SE |\n|---:|---|---:|---:|")
        for r in block_fits
            @printf(io,"| %.0f | %s | %d | %.6f ± %.6f |\n",r["q"],r["observable"],r["blocksize"],r["exponent"],r["stderr"])
        end
        println(io)
        shifts=Float64[]
        for row in block_fits
            baseline=only(filter(r->r["q"]==row["q"] && r["observable"]==row["observable"] && r["blocksize"]==config["primary_block"],block_fits))
            push!(shifts,abs(row["exponent"]-baseline["exponent"])/baseline["stderr"])
        end
        @printf(io,"The largest change of fitted slope across these block choices is %.2f times the primary fit SE. This limited block sensitivity does not explain away the individual-chain flags; additional long independent runs would make the uncertainty estimates more stable.\n\n",maximum(shifts))
        println(io,"Three restart audits (including L=64 bond and spin runs) reproduce observations, final configurations, work/clock counters and future RNG. Metadata contains seeds, units, Julia version and source/input hashes. Reanalysis never resamples physical states. IBFS improves connectivity searches but is not a fully dynamic connectivity data structure; no dynamic critical exponent or sub-sweep acceleration is inferred.\n")
        for (name,title) in (("critical_scaling","Scaling and scaled amplitudes"),("fit_diagnostics","Cutoff sensitivity and effective exponents"),("connectivity_performance","Paired connectivity benchmark"),("chain_diagnostics","Correlations and representation checks"))
            println(io,"\n![$title](figures/$name.png)")
        end
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    analyze(isempty(ARGS) ? joinpath(@__DIR__,"results","critical") : abspath(ARGS[1]))
end
