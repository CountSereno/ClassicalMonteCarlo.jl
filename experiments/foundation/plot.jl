# Julia-native offline figures. Run with the separate plotting environment.
using CairoMakie, Statistics, TOML, SHA
include("tables.jl")
using .FoundationTables

const NAMES = Dict("metropolis"=>"Metropolis","wolff"=>"Wolff","swendsen_wang"=>"Swendsen–Wang")
const SHORT_NAMES = Dict("metropolis"=>"Metro.","wolff"=>"Wolff","swendsen_wang"=>"SW")
const COLORS = Dict("metropolis"=>"#2864a0","wolff"=>"#d07b35","swendsen_wang"=>"#248b81")
number(row,key) = parse(Float64,row[key])
selected(rows;conditions...) = filter(r->all(r[string(k)]==string(v) for (k,v) in conditions),rows)

function save_figure(dir,name,figure)
    save(joinpath(dir,name*".png"),figure;px_per_unit=1.8)
    save(joinpath(dir,name*".svg"),figure)
    save(joinpath(dir,name*".pdf"),figure)
end

function plot_all(root)
    tables = joinpath(root,"analysis")
    aggregate=read_table(joinpath(tables,"aggregate.csv"))
    chains=read_table(joinpath(tables,"chains.csv"))
    blocks=read_table(joinpath(tables,"block_errors.csv"))
    correlations=read_table(joinpath(tables,"acf.csv"))
    cfg=TOML.parsefile(joinpath(root,"config.toml"))
    algorithms=cfg["algorithms"]
    nchains=length(cfg["initial_states"])
    figdir=joinpath(root,"figures");mkpath(figdir)
    set_theme!(Theme(fontsize=14,Axis=(xgridvisible=false,ygridcolor=(:black,.07),
        topspinevisible=false,rightspinevisible=false),Legend=(framevisible=false,)))

    fig=Figure(size=(1200,720))
    Label(fig[0,1:3],"Small-system validation against independent enumeration",fontsize=23,font=:bold)
    Label(fig[3,1:3],"Periodic L=$(cfg["exact_L"]), T=$(cfg["exact_T"]); $nchains chains per algorithm · bars: diagnostic ±1 SE · dashed: exact",fontsize=13,color=:gray35)
    observables=[("energy_density","Energy density","E / N"),
        ("abs_magnetization_density","Absolute magnetization","⟨|M|⟩ / N"),
        ("heat_capacity","Heat capacity per site","cᵥ"),
        ("susceptibility","Signed-M susceptibility","χ"),
        ("susceptibility_abs","Absolute-M fluctuation proxy","χ_abs"),
        ("binder","Binder cumulant","U₄")]
    for (index,(key,title,ylabel)) in enumerate(observables)
        axis=Axis(fig[(index-1)÷3+1,(index-1)%3+1];title,ylabel,
            xticks=(1:length(algorithms),[SHORT_NAMES[a] for a in algorithms]),xticklabelrotation=0)
        entries=selected(aggregate;experiment="exact",observable=key)
        reference=number(first(entries),"exact")
        hlines!(axis,[reference];color=:gray25,linestyle=:dash,linewidth=1.4)
        for (x,algorithm) in enumerate(algorithms)
            row=only(selected(entries;algorithm))
            value=number(row,"mean");se=number(row,"diagnostic_se")
            errorbars!(axis,[x],[value],[se];color=COLORS[algorithm],whiskerwidth=12,linewidth=2)
            scatter!(axis,[x],[value];color=COLORS[algorithm],markersize=11)
        end
        xlims!(axis,.5,length(algorithms)+.5)
    end
    save_figure(figdir,"exact_validation",fig)

    fig=Figure(size=(1300,490))
    Label(fig[0,1:3],"Critical Ising: observable correlation and measured implementation cost",fontsize=22,font=:bold)
    axes=[Axis(fig[1,i];xlabel="Linear size L",yscale=log10,xticks=cfg["sizes"],
        title=t,ylabel=y) for (i,(t,y)) in enumerate([
        ("Native algorithm clock","Energy τ in algorithm updates"),
        ("Approximate spin-work clock","Energy τ in spin sweeps"),
        ("Includes production measurements","Estimated energy ESS / second")])]
    keys=["tau_energy_updates","tau_energy_spin_sweeps","ESS_energy_per_second"]
    for algorithm in algorithms
        for (axis,key) in zip(axes,keys)
            rows=selected(chains;experiment="critical",algorithm)
            x=cfg["sizes"]
            y=[mean(number(r,key) for r in selected(rows;L)) for L in x]
            lines!(axis,x,y;color=COLORS[algorithm],linewidth=2,label=NAMES[algorithm])
            scatter!(axis,x,y;color=COLORS[algorithm],markersize=9)
            for row in rows
                scatter!(axis,[number(row,"L")],[number(row,key)];color=(COLORS[algorithm],.3),markersize=5)
            end
        end
    end
    Legend(fig[2,1:3],axes[1];orientation=:horizontal)
    Label(fig[3,1:3],"Lines: $nchains-chain means; dots: individual chains · Wolff spin-work scaling uses mean cluster size / N · no exponent fit",fontsize=12,color=:gray35)
    truncated=count(row->row["experiment"]=="critical" && parse(Int,row["window_energy"])==cfg["maxlag"],chains)
    Label(fig[4,1:3],"Primary maximum lag: $(cfg["maxlag"]); $truncated energy chains reach this limit · τ may be low and ESS high; see window-sensitivity report",fontsize=12,color=:gray35)
    save_figure(figdir,"algorithm_efficiency",fig)

    fig=Figure(size=(1180,780))
    Lmax=maximum(cfg["sizes"])
    Label(fig[0,1:2],"Correlation and block-length diagnostics at L=$Lmax",fontsize=23,font=:bold)
    axE=Axis(fig[1,1];title="Energy autocorrelation",xlabel="Lag in stored algorithm updates",ylabel="ρ_E(t)")
    axM=Axis(fig[1,2];title="Absolute-magnetization autocorrelation",xlabel="Lag in stored algorithm updates",ylabel="ρ_|M|(t)")
    axBE=Axis(fig[2,1];title="Energy mean: block jackknife error",xlabel="Block size (stored records)",ylabel="Within-chain standard error",xscale=log2,yscale=log10,xticks=cfg["block_sizes"])
    axBU=Axis(fig[2,2];title="Binder: paired-block jackknife error",xlabel="Block size (stored records)",ylabel="Within-chain standard error",xscale=log2,yscale=log10,xticks=cfg["block_sizes"])
    for algorithm in algorithms
        for (axis,observable) in ((axE,"energy"),(axM,"absM"))
            rows=selected(correlations;algorithm,observable)
            lags=sort(unique(number(r,"lag") for r in rows))
            means=[mean(number(r,"acf") for r in rows if number(r,"lag")==lag) for lag in lags]
            for chain in 1:length(cfg["initial_states"])
                cr=selected(rows;chain)
                lines!(axis,[number(r,"lag") for r in cr],[number(r,"acf") for r in cr];color=(COLORS[algorithm],.2),linewidth=.8)
            end
            lines!(axis,lags,means;color=COLORS[algorithm],linewidth=2,label=NAMES[algorithm])
        end
        for (axis,observable) in ((axBE,"energy_density"),(axBU,"binder"))
            rows=selected(blocks;experiment="critical",L=Lmax,algorithm,observable)
            bs=cfg["block_sizes"]
            for chain in 1:length(cfg["initial_states"])
                cr=selected(rows;chain)
                vals=[number(only(selected(cr;blocksize=b)),"stderr") for b in bs]
                lines!(axis,bs,vals;color=(COLORS[algorithm],.2),linewidth=.8)
            end
            vals=[mean(number(r,"stderr") for r in selected(rows;blocksize=b)) for b in bs]
            lines!(axis,bs,vals;color=COLORS[algorithm],linewidth=2)
            scatter!(axis,bs,vals;color=COLORS[algorithm],markersize=8)
        end
    end
    for axis in (axE,axM)
        hlines!(axis,[0.];color=:gray50,linewidth=.8)
        xlims!(axis,0,min(250,cfg["maxlag"]))
        ylims!(axis,-.1,1.02)
    end
    Legend(fig[3,1:2],axE;orientation=:horizontal)
    Label(fig[4,1:2],"Faint curves: individual chains; bold curves: means · a flat block-error curve is a diagnostic, not proof of equilibration",fontsize=12,color=:gray35)
    save_figure(figdir,"correlation_diagnostics",fig)
    plotenv=normpath(joinpath(@__DIR__,"../plotting"))
    info=Dict("julia_version"=>string(VERSION),"plot_source_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "table_reader_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"tables.jl")))),
        "cairomakie_version"=>string(Base.pkgversion(CairoMakie)),
        "input_sha256"=>Dict(f=>bytes2hex(sha256(read(joinpath(tables,f)))) for f in ("aggregate.csv","chains.csv","block_errors.csv","acf.csv")))
    if isfile(joinpath(plotenv,"Manifest.toml"))
        info["plot_manifest_sha256"]=bytes2hex(sha256(read(joinpath(plotenv,"Manifest.toml"))))
    end
    info["plot_project_sha256"]=bytes2hex(sha256(read(joinpath(plotenv,"Project.toml"))))
    open(io->TOML.print(io,info),joinpath(figdir,"metadata.toml"),"w")
    println("Wrote PNG, SVG and PDF figures to $figdir")
end

if abspath(PROGRAM_FILE)==@__FILE__
    plot_all(length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(@__DIR__,"results","demo"))
end
