# 模型、概率和时间单位

## 自旋模型

二维矩形周期方格，每个边长至少为 3。每个正方向最近邻键计算一次，令 k_B=1、β=1/T，无外场：

- Ising：s_i ∈ {−1,+1}，H=−J Σ⟨ij⟩ s_i s_j。`energy` 返回总能量，`magnetization` 返回总有符号磁化 M=Σs_i。
- Potts：s_i ∈ {1,…,q}，H=−J Σ⟨ij⟩ δ(s_i,s_j)。`magnetization` 返回 (q max_a N_a/N−1)/(q−1)，范围 0 到 1，不是总磁化。二态 Potts 的 J 与 Ising 的 J 有因子 2 的换算。

Metropolis 支持任意实数 J。一个 sweep 是 N 次随机选点、有放回的提议，不是顺序遍历。`metropolis!` 返回接受次数。Wolff 和 SW 要求 J≥0：对相同自旋的键，Ising 和 Potts 激活概率分别为 1−exp(−2βJ) 和 1−exp(−βJ)。Wolff 返回更新簇大小，一次调用是一次簇更新，不是 sweep；SW 重新着色所有簇，返回簇数。在 βJ=0 的 Ising 或二态 Potts 情况，为避免退化周期性，Metropolis 对每次提议、Wolff 对每次调用均有 1/2 概率空操作；此时 Wolff 可能返回 0。时间尺度不能直接比较而不说明归一化。有限精度下 βJ 极大时激活概率可能舍入到 1，更新可能丧失遍历性。

## 键渗流与 Sweeny

`BondGraph` 是无重边、无自环的有限无向简单图。开放方格支持 L≥1；周期方格 L≥3。独立键渗流对每条键独立以概率 p 占据；`cluster_sizes` 包括孤立点。当前没有 site percolation、跨越或绕环观测量。

Sweeny 随机团簇分布为 `p^b (1-p)^(m-b) q^k`，b 是占据键数，m 是总键数，k 包含孤立点在内的连通分量数，0<p<1、q>0。每次从所有键均匀选一条，以 1/2 概率保持不变，否则对翻转做 Metropolis 接受/拒绝；一个 sweep 做 m 次这样的有放回尝试。每次更新重建并查集，删除桥的判断正确但速度不适合大系统。

## Worm 与 Wang–Landau

`IsingWorm` 对零场铁磁 Ising 使用高温图展开，x=tanh(βJ)，图权重为 x^b。尾端固定，头端每步均匀选邻点，提议翻转该键；添加接受率 x，删除接受率 1。周期 L≥3，且要求在浮点表示中 0<x<1。当前状态是图而不是自旋构型；占据图除头尾外所有顶点度数为偶数。

按每一次**尝试**测量（包括拒绝），闭合比例的倒数估计零场无量纲量 ⟨M²⟩/N；物理磁化率按通常约定为 β⟨M²⟩/N。为误差分析保存闭合指示序列并使用适当分块或独立链；倒数是非线性统计量。不要仅以闭合事件计时后套用这一估计式。

Wang–Landau 的 `bin` 需将所有可达能量或其他宏观变量映射到连续的 1:nbins；不可达的 bin 会阻止平坦性判断。`propose` 必须对称、返回独立的新状态且不能原地更改输入；边界处改变提议概率会破坏对称性。输出的 log 密度只确定到一个加法常数。直方图达到平坦标准并完成修改因子递减只代表该停止规则被满足，并不保证物理精度。示例是只有三个微观状态的可解析玩具模型，真实能量范围、遍历性和误差需要单独验证。

## 误差和有限尺寸

分析工具假设输入已平衡且平稳。分块要求长度整除，且块长要足够大到使块近似独立。`jackknife` 是删除连续块的标量统计量 jackknife。自相关使用简单的首个非正值截断；噪声较大或强负相关的链需要更审慎的窗口选择。`finite_size_fit` 仅拟合预先给定 ω 的 `y(L)=a+b L^(-ω)`；加权拟合假定所给标准差已知且样本点独立。它不是通用的临界点搜索、数据坍缩或全套 FSS 分析。

`ising_observables(E, M; beta, nsites, blocksize)` 对配对的**总量**序列同时删除相同的时间块，保留 E、M 及高阶矩的协方差。令 N 为格点数，各个 `raw` 定义为：

- 能量密度 ⟨E⟩/N，绝对磁化密度 ⟨|M|⟩/N；
- 每格点比热 C/N = β²(⟨E²⟩−⟨E⟩²)/N；
- 每格点磁化率 χ/N = β(⟨M²⟩−⟨M⟩²)/N；零场下仍需充分遍历两个磁化扇区；
- `susceptibility_abs` = β(⟨M²⟩−⟨|M|⟩²)/N，仅是另一种涨落量，不能混称有限系统的线性响应磁化率；
- U₄ = 1−⟨M⁴⟩/(3⟨M²⟩²)。

各字段提供 `raw`、偏差修正后的 `estimate`、`bias` 和 `stderr`。偏差修正量在短样本下可能超出物理区间；块长和热化需要独立诊断。新增实验保留多个块长、独立初态、原始时间序列和警告，见 [实验说明](../experiments/README.md)。

## Potts 与随机团簇的配对统计

`potts_observables(E, order; beta, nsites, blocksize, simplex_squared=nothing)`
使用总能量 E 与无量纲最大颜色序参量
`order=(q max_a N_a/N−1)/(q−1)`。返回的定义是：

- `energy_density` = ⟨E⟩/N；`order_parameter` = ⟨order⟩；
- `heat_capacity` = β² Var(E)/N；
- `order_parameter_fluctuation` = βN Var(order)，是最大颜色序参量的涨落代理量，不能当作某个固定磁场方向的线性响应磁化率；
- `moment_ratio` = ⟨order⁴⟩/⟨order²⟩²，不套用 Ising 标量 Binder 累积量的系数；若完整样本或任一删块样本的 ⟨order²⟩ 为零，该比值未定义，函数明确报错。

可选的 `simplex_squared` 是 Potts simplex 向量序参量的平方：
`(q Σ_a N_a²/N²−1)/(q−1)`，其平均值作为 `simplex_second_moment` 返回。
对于 q>2，它与最大颜色序参量的平方不同。该序列与 E、order 使用完全相同的时间块，
不会被单独打乱或重采样。能量和序参量方差先中心化再计算，减少能量零点较大时的消减误差。

`random_cluster_observables(b,k,cmax,s2; nedges,nvertices,blocksize)`
接受同一次测量的占据键数 b、连通分量数 k、最大分量大小 cmax，以及
`s2=Σ_C |C|²`。所有孤立点与最大团簇都计入。记总边数为 m、顶点数为 N，则：

- `bond_density` = ⟨b⟩/m；`cluster_density` = ⟨k⟩/N；
- `largest_cluster_fraction` = ⟨cmax⟩/N；
- `mean_cluster_size` = ⟨s2⟩/N，是均匀选取一个顶点后，它所属团簇的平均大小，包含最大团簇；
- `bond_variance_per_edge` = Var(b)/m，是几何涨落量，未自动解释为比热。

该统计接口要求 m>0，因为无边图的键密度没有定义。p=0 或 p=1 下的常量观测序列可以分析，
给出零样本方差；零方差不代表任意系综下都已充分采样。以上接口返回与 Ising 统计一致的
`raw/estimate/bias/stderr` 和分块元数据，使用共同的配对矩删块计算，保持非线性量误差中的协方差。

零场铁磁整数 q Potts 与 FK 随机团簇可用 `p=1−exp(−βJ)` 联系。在相同有限图上，p>0 时有
`⟨E_Potts⟩/N = −J⟨b⟩/(pN)`，以及
`⟨simplex_squared⟩ = ⟨s2⟩/N²`。第二个关系来自逐簇独立均匀着色，
是对系综平均的等式，不能要求两条独立轨迹逐个构型相等。非整数 q 的随机团簇本身是合法概率模型，
但不能因此声称它具有整数颜色自旋表示。p=0 时能量关系不能直接除以 p。

[多模型实验](../experiments/multimodel/README.md) 用独立颜色枚举与 BFS 键子集枚举检查这些关系；
精确枚举不调用生产代码中的能量或并查集函数。
