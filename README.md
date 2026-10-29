# ClassicalMonteCarlo.jl

**Reproducible Monte Carlo experiments for classical statistical physics.**
This Julia package implements Ising/Potts spin updates, random-cluster sampling,
worm sampling, and basic uncertainty analysis. Small-system exact checks and
explicit model conventions take priority over feature count. It is an educational
research-software foundation; the included short runs do not establish precision
critical exponents.

- [Numerical experiments, raw traces and figures](experiments/README.md)
- [First development cycle: exact checks, correlation diagnostics and Julia figures](experiments/foundation/README.md)
- [Software design](docs/software-design.md) · [Pause, resume and export](docs/runner.md)
- [Model and statistical conventions](docs/model-conventions.md)
- [Tests](test/runtests.jl) · [Development roadmap](docs/roadmap.md)

一个从零构建的 Julia 经典统计物理 Monte Carlo 工具包。当前版本以可读、可复现和小系统可验证为优先，既可用来学习算法，也可作为扩展研究代码的基础。它仍是早期版本，不应未经独立验证便用于生产规模计算或发表结果。

## 已实现

| 部件 | 内容与范围 |
| --- | --- |
| 自旋模型 | 二维周期方格 Ising 和整数 q 态 Potts，零外场 |
| 自旋更新 | Metropolis、单簇 Wolff、Swendsen–Wang；后两者限铁磁耦合 |
| 渗流与随机团簇 | 独立键渗流、一般有限简单图、组件大小、非整数 q 的 Sweeny 更新 |
| Worm | 二维周期方格 Ising 高温图展开，固定尾端的局域 worm |
| Entropic sampling | 通用 Wang–Landau 平直直方图采样器，调用者提供状态、能量分箱和对称提议 |
| 工具 | 独立 RNG、连续分块、联合热力学量块 jackknife、自相关时间和固定修正指数的有限尺寸外推 |
| 可移植核心 | 独立 Metropolis–Hastings 接受过滤器；Wolff/SW 可复用工作区 |
| 实验运行 | Ising 专用的参数校验、显式更新时钟、暂停恢复、CSV/TOML 导出和源码指纹 |
| 质量保障 | 小系统精确枚举、随机团簇细致平衡、边界条件、可复现性测试，以及 GitHub Actions 配置 |

核心包只使用 Julia 标准库，支持 Julia 1.10 及以上；本机实际验证版本见下。CairoMakie 绘图依赖保存在独立的实验环境中。各功能的定义、限制和时间单位见 [模型与算法约定](docs/model-conventions.md)。参考文献和本机资料见 [参考资料](docs/references.md)，进一步的开发路径见 [路线图](docs/roadmap.md)。

## 快速开始

克隆或下载仓库后，在仓库根目录执行：

```bash
julia --project=. examples/ising.jl
julia --project=. examples/random_cluster.jl
julia --project=. examples/worm.jl
julia --project=. examples/wang_landau.jl
julia --project=. examples/acceptance.jl
```

例如：

```julia
using ClassicalMonteCarlo, Statistics

rng = new_rng(20261028)
m = Ising(16; init=:random, rng=rng)
beta = log(1 + sqrt(2)) / 2 # 方格 Ising 的临界逆温，J=1
for _ in 1:1000
    wolff!(m, beta; rng=rng) # 热化：一次调用更新一个簇
end
samples = Float64[]
for _ in 1:10000
    wolff!(m, beta; rng=rng)
    push!(samples, energy(m) / length(m.spins))
end
result = integrated_autocorrelation_time(samples; maxlag=500)
println((mean=mean(samples), stderr=result.stderr_mean, tau=result.tau))
```

反复使用簇更新时，可先构造 `workspace = WolffWorkspace(m)`，再传入
`wolff!(m, beta; rng=rng, workspace=workspace)` 复用内存。需要完整的运行记录与恢复时：

```julia
spec = RunSpec(L=16, beta=beta, algorithm=:wolff, seed=20261030,
               burn=2000, samples=8000, stride=1)
state = start_run(spec)
advance!(state, 2500)
save_checkpoint("results/example.checkpoint", state)
state = load_checkpoint("results/example.checkpoint"; expected_spec=spec)
run!(state)
export_run("results/example", state) # 新目录，保存配对 E/M 序列与元数据
```

[本轮实验](experiments/foundation/README.md) 将生成、验证、统计分析和绘图分成独立命令。
其精确对照与临界点算法比较用于检验实现和分析流程，不用于精确临界指数估计。

样例运行参数仅供演示。要开展可靠计算，请检查热化、不同初态与独立随机种子的结果、自相关时间与分块稳定性，再选择模拟长度和误差估计。单次短链或一个收敛标志不能证明平衡。记录种子、Julia 版本、代码版本、模型参数及更新和测量间隔；Julia 的 RNG 流不保证跨版本逐位一致。

## 测试

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

没有 Julia package registry 或不希望让 Pkg 修改全局 depot 时，标准库依赖可直接运行测试：

```bash
julia --project=. test/runtests.jl
```

本次在 Julia 1.12.6 上执行完整测试及示例。CI 也配置了 Julia 1.10 和当前稳定版，但需要 GitHub Actions 实际运行后才能确认对应平台的结果。

## 接口一览

- `new_rng(seed)`、`independent_rngs(seed, n)`；返回独立使用的标准库 Xoshiro 生成器。
- `metropolis_accept(logratio, u)`、`metropolis_accept(rng, logratio)`、`metropolis_logratio(beta, delta_energy; log_proposal_ratio)`；可脱离其他模块独立使用。
- `Ising(Lx, Ly=Lx; ...)`、`Potts(Lx, Ly, q; ...)`；`energy`、`magnetization`。
- `metropolis!(model, beta; rng, sweeps)`、`wolff!(model, beta; rng)`、`swendsen_wang!(model, beta; rng)`。
- `WolffWorkspace(model)`、`SwendsenWangWorkspace(model)`；通过更新函数的 `workspace` 关键字复用，每条链独占。
- `RunSpec`、`start_run`、`advance!`、`run!`、`save_checkpoint`、`load_checkpoint`、`export_run`；当前运行器限方格 Ising。
- `square_bond_graph`、`BondConfiguration`、`sample_bond_percolation!`、`cluster_sizes`、`nclusters`、`sweeny_update!`、`sweeny_sweep!`。
- `IsingWorm`、`worm_step!`、`worm_closed`、`worm_bond_count`。
- `wang_landau(rng, initial; nbins, bin, propose, ...)`。
- `block_means`、`jackknife`、`autocorrelation`、`integrated_autocorrelation_time`、`finite_size_fit`。
- `ising_observables(E, M; beta, nsites, blocksize)`：输入总能量与有符号总磁化的配对序列，返回能量密度、绝对磁化密度、比热、两种磁化涨落量和 Binder 累积量的联合分块 jackknife 估计。

使用 Julia 内置帮助，如 `?wang_landau`，查看参数和详细假设。源码和测试是约定的一部分。

## 许可证

本项目采用 [MIT 许可证](LICENSE)，版权署名为 Dieucra。
