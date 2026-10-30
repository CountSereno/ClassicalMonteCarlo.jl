# Run with the separately pinned experiments/plotting environment.
using CairoMakie, Statistics, TOML, SHA
include("common.jl")
include("../multimodel/tables.jl")
using .MultiModelTables

number(row,key)=parse(Float64,row[key])
column(rows,key)=[number(row,key) for row in rows]
choose(rows;kwargs...)=filter(row->all(row[string(k)]==string(v) for (k,v) in kwargs),rows)
ordered(rows,key)=sort(rows;by=row->number(row,key))
const COLOURS=Dict(2.0=>"#27649A",3.0=>"#CB6032")
const SIZE_TICKS=([8,16,32,64],["8","16","32","64"])
const OBS_LABELS=Dict("C1"=>"Largest FK cluster C₁","S2_over_N"=>"Cluster second moment S₂ / N")
reference_power(q,obs)=obs=="C1" ? (q==2 ? 15/8 : 28/15) : (q==2 ? 7/4 : 26/15)

function save_figure(fig,root,name)
    for extension in ("png","svg","pdf")
        path=joinpath(root,"figures",name*"."*extension)
        extension=="png" ? save(path,fig;px_per_unit=1.5) : save(path,fig)
    end
end

function plot_scaling(root,aggregate,fits)
    fig=Figure(size=(1240,920),fontsize=18)
    Label(fig[0,1:2],"Critical FK clusters: seven sizes and independent chains",fontsize=25,font=:bold)
    for (col,obs) in enumerate(("C1","S2_over_N"))
        raw=Axis(fig[1,col],title=OBS_LABELS[obs],xlabel="Linear size L",ylabel="Mean ± diagnostic SE",
            xscale=log2,yscale=log10,xticks=SIZE_TICKS)
        scaled=Axis(fig[2,col],xlabel="Linear size L",ylabel="Mean / L^(reference power)",xscale=log2,xticks=SIZE_TICKS)
        for q in (2.,3.)
            rows=ordered(choose(aggregate;family="RandomCluster",q,observable=obs),"L")
            x,y,se=column(rows,"L"),column(rows,"estimate"),column(rows,"stderr")
            colour=COLOURS[q]
            errorbars!(raw,x,y,se;color=colour,whiskerwidth=8)
            scatterlines!(raw,x,y;color=colour,label="q=$(Int(q))")
            fit=only(choose(fits;q,observable=obs,Lmin=16.0))
            grid=range(16,64;length=150)
            lines!(raw,grid,number(fit,"amplitude").*grid.^number(fit,"exponent");color=colour,linestyle=:dash)
            scale=x.^reference_power(q,obs)
            errorbars!(scaled,x,y./scale,se./scale;color=colour,whiskerwidth=8)
            scatterlines!(scaled,x,y./scale;color=colour,label="q=$(Int(q))")
        end
        axislegend(raw;position=:lt)
    end
    Label(fig[3,1:2],"Periodic square · p = √q / (1 + √q) · dashed: leading-power fit with L ≥ 16",fontsize=16,color=:gray40)
    Label(fig[4,1:2],"Reference powers: C₁ = 15/8, 28/15; S₂/N = 7/4, 26/15 · corrections and cutoff sensitivity remain",fontsize=16,color=:gray40)
    save_figure(fig,root,"critical_scaling")
end

function plot_fits(root,fits,effective)
    fig=Figure(size=(1240,920),fontsize=18)
    Label(fig[0,1:2],"Finite-size sensitivity: statistical errors do not include fit-model bias",fontsize=24,font=:bold)
    for (col,obs) in enumerate(("C1","S2_over_N"))
        cutoff=Axis(fig[1,col],title=OBS_LABELS[obs],xlabel="Smallest included size Lmin",ylabel="Leading-power fitted exponent",xticks=([8,12,16,24,32],["8","12","16","24","32"]))
        pairs=Axis(fig[2,col],xlabel="Geometric mean of adjacent sizes",ylabel="Adjacent-size effective exponent",xscale=log2,xticks=([10,16,24,32,48],["10","16","24","32","48"]))
        for q in (2.,3.)
            colour=COLOURS[q]
            rows=ordered(choose(fits;q,observable=obs),"Lmin")
            errorbars!(cutoff,column(rows,"Lmin"),column(rows,"exponent"),column(rows,"stderr");color=colour,whiskerwidth=8)
            scatterlines!(cutoff,column(rows,"Lmin"),column(rows,"exponent");color=colour,label="q=$(Int(q))")
            hlines!(cutoff,[reference_power(q,obs)];color=(colour,.65),linestyle=:dash)
            rows=ordered(choose(effective;q,observable=obs),"L_geometric")
            errorbars!(pairs,column(rows,"L_geometric"),column(rows,"exponent"),column(rows,"stderr");color=colour,whiskerwidth=8)
            scatterlines!(pairs,column(rows,"L_geometric"),column(rows,"exponent");color=colour)
            hlines!(pairs,[reference_power(q,obs)];color=(colour,.65),linestyle=:dash)
        end
        axislegend(cutoff;position=:lt)
    end
    Label(fig[3,1:2],"Dashed lines: theoretical powers · all fits end at L=64 · Lmin=32 has only one residual degree of freedom",fontsize=16,color=:gray40)
    Label(fig[4,1:2],"Adjacent estimates share one size and are anticorrelated; their full covariance is exported",fontsize=16,color=:gray40)
    save_figure(fig,root,"fit_diagnostics")
end

function plot_performance(root,benchmark)
    fig=Figure(size=(1240,600),fontsize=18)
    Label(fig[0,1:2],"Sweeny connectivity: paired trajectories, different search cost",fontsize=25,font=:bold)
    timing=Axis(fig[1,1],xlabel="Linear size L",ylabel="Milliseconds per sweep",xscale=log2,yscale=log10,xticks=SIZE_TICKS)
    speed=Axis(fig[1,2],xlabel="Linear size L",ylabel="Reference / IBFS elapsed time",xscale=log2,xticks=SIZE_TICKS)
    for q in (2.,3.)
        rows=ordered(choose(benchmark;q),"L"); x=column(rows,"L"); colour=COLOURS[q]
        scatterlines!(timing,x,1000column(rows,"reference_seconds");color=colour,linestyle=:dash,label="q=$(Int(q)) rebuild")
        scatterlines!(timing,x,1000column(rows,"ibfs_seconds");color=colour,label="q=$(Int(q)) IBFS")
        scatterlines!(speed,x,column(rows,"speedup");color=colour,label="q=$(Int(q))")
    end
    axislegend(timing;position=:lt); axislegend(speed;position=:lt)
    Label(fig[2,1:2],"Median of five alternating-order repeats · update only · same configurations and future RNG in each pair",fontsize=16,color=:gray40)
    Label(fig[3,1:2],"Both reuse scratch and skip lazy queries · 0 allocated bytes per sweep after warmup · single-machine timings",fontsize=16,color=:gray40)
    save_figure(fig,root,"connectivity_performance")
end

function plot_diagnostics(root,acf,block_errors,cross)
    fig=Figure(size=(1240,920),fontsize=18)
    Label(fig[0,1:2],"Correlation, block sensitivity and representation checks",fontsize=25,font=:bold)
    a=Axis(fig[1,1],title="L=64 second-moment autocorrelation",xlabel="Lag in recorded sweeps",ylabel="ACF")
    b=Axis(fig[1,2],title="All Sweeny chains: second-moment block errors",xlabel="Linear size L",ylabel="2048-block SE / 512-block SE",xscale=log2,
        xticks=SIZE_TICKS)
    for q in (2.,3.)
        colour=COLOURS[q]
        rows=choose(acf;family="RandomCluster",q,L=64,observable="S2_over_N")
        ids=sort!(unique(row["run_id"] for row in rows))
        series=Vector{Float64}[]
        for id in ids
            chain=ordered(filter(row->row["run_id"]==id && number(row,"lag")<=300,rows),"lag")
            values=column(chain,"acf"); push!(series,values)
            lines!(a,column(chain,"lag"),values;color=(colour,.18))

        end
        lines!(a,0:length(first(series))-1,mean(series);color=colour,linewidth=3,label="q=$(Int(q)) Sweeny")
        all_blocks=choose(block_errors;family="RandomCluster",q,observable="S2_over_N")
        sizes=sort!(unique(number(row,"L") for row in all_blocks))
        averages=Float64[]
        for L in sizes
            group=filter(row->number(row,"L")==L,all_blocks)
            ratios=Float64[]
            for id in unique(row["run_id"] for row in group)
                numerator=number(only(choose(group;run_id=id,blocksize=2048)),"stderr")
                denominator=number(only(choose(group;run_id=id,blocksize=512)),"stderr")
                push!(ratios,numerator/denominator)
            end
            scatter!(b,fill(L,length(ratios)),ratios;color=(colour,.4))
            push!(averages,mean(ratios))
        end
        lines!(b,sizes,averages;color=colour,linewidth=2)
    end
    hlines!(b,[1.];color=:gray50,linestyle=:dash); hlines!(b,[1.5];color=:gray50,linestyle=:dot); axislegend(a;position=:rt)
    for (col,obs) in enumerate(("energy_density","S2_over_N"))
        ax=Axis(fig[2,col],title=obs=="energy_density" ? "q=3 energy identity" : "q=3 simplex / FK identity",
            xlabel="Linear size L",ylabel="(Spin − bond) / combined SE",xticks=([16,32,64],["16","32","64"]))
        rows=ordered(choose(cross;observable=obs),"L")
        hspan!(ax,-1,1;color=(:gray,.12)); hlines!(ax,[0.];color=:gray50,linestyle=:dash)
        scatterlines!(ax,column(rows,"L"),column(rows,"z");color="#7354A3",markersize=12)
    end
    Label(fig[3,1:2],"Faint curves: four independent chains · block checks and finite-window ACF do not establish equilibrium",fontsize=16,color=:gray40)
    Label(fig[4,1:2],"Spin SW and bond Sweeny use independent seeds · shaded band: ±1 combined diagnostic SE",fontsize=16,color=:gray40)
    save_figure(fig,root,"chain_diagnostics")
end

function plot_study(root)
    set_theme!(Theme(Axis=(xgridcolor=(:gray,.13),ygridcolor=(:gray,.13),)))
    mkpath(joinpath(root,"figures"))
    table(name)=read_table(joinpath(root,"analysis",name*".csv"))
    plot_scaling(root,table("aggregate"),table("power_fits"))
    plot_fits(root,table("power_fits"),table("effective_exponents"))
    benchmark_path=joinpath(REPOSITORY,"benchmarks/results/connectivity/summary.csv")
    plot_performance(root,read_table(benchmark_path))
    plot_diagnostics(root,table("acf"),table("block_errors"),table("fk_comparison"))
    inputs=["aggregate","power_fits","effective_exponents","acf","block_errors","fk_comparison"]
    metadata=Dict("julia_version"=>string(VERSION),"cairomakie_version"=>string(Base.pkgversion(CairoMakie)),
        "source_sha256"=>source_hashes(["experiments/scaling/plot.jl","experiments/scaling/common.jl",
            "experiments/multimodel/tables.jl","experiments/plotting/Project.toml","experiments/plotting/Manifest.toml"]),
        "input_sha256"=>Dict("analysis/"*name*".csv"=>file_hash(joinpath(root,"analysis",name*".csv")) for name in inputs),
        "benchmark_sha256"=>file_hash(benchmark_path),
        "output_sha256"=>Dict(name=>file_hash(joinpath(root,"figures",name)) for name in readdir(joinpath(root,"figures")) if splitext(name)[2] in (".png",".svg",".pdf")))
    save_toml(joinpath(root,"figures","metadata.toml"),metadata)
end

if abspath(PROGRAM_FILE)==@__FILE__
    plot_study(isempty(ARGS) ? joinpath(@__DIR__,"results","critical") : abspath(ARGS[1]))
end
