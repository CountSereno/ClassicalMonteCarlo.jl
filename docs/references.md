# 参考资料与来源

实现根据公开的算法定义独立编写；参考资料不意味着作者认可本项目。论文、示例代码分别受其各自许可证及版权条款约束。

## 经典算法与验证

- Metropolis et al., “Equation of State Calculations by Fast Computing Machines” (1953), [DOI: 10.1063/1.1699114](https://doi.org/10.1063/1.1699114).
- Swendsen and Wang, “Nonuniversal critical dynamics in Monte Carlo simulations” (1987), [DOI: 10.1103/PhysRevLett.58.86](https://doi.org/10.1103/PhysRevLett.58.86).
- Wolff, “Collective Monte Carlo Updating for Spin Systems” (1989), [DOI: 10.1103/PhysRevLett.62.361](https://doi.org/10.1103/PhysRevLett.62.361).
- Sweeny, “Monte Carlo study of weighted percolation clusters relevant to the Potts models” (1983), [DOI: 10.1103/PhysRevB.27.4445](https://doi.org/10.1103/PhysRevB.27.4445).
- Beffara and Duminil-Copin, “The self-dual point of the two-dimensional random-cluster model is critical for q≥1”, [arXiv:1006.5073](https://arxiv.org/abs/1006.5073). 多模型实验代入方格的已知临界参数，未重新拟合临界点。
- Elçi and Weigel, “Efficient simulation of the random-cluster model” (2013), [arXiv:1307.6647](https://arxiv.org/abs/1307.6647). 参考 FK 表示、几何观测量和连通性成本的区分；本包尚未实现该工作的动态连通性算法。
- Deng, Garoni and Sokal, “Critical speeding-up in a local dynamics for the random-cluster model” (2007), [arXiv:cond-mat/0701113](https://arxiv.org/abs/cond-mat/0701113). 本轮每 sweep 测量一次，不能分辨其研究中的子 sweep 快过程，亦未测定动力学指数。
- Prokof’ev and Svistunov, “Worm Algorithms for Classical Statistical Models” (2001), [arXiv](https://arxiv.org/abs/cond-mat/0103146).
- Wang and Landau, “Determining the density of states for classical statistical models” (2001), [arXiv](https://arxiv.org/abs/cond-mat/0107006).
- Newman and Ziff, [percolation algorithm and original implementation](https://websites.umich.edu/~mejn/percolation/). A useful future implementation and validation target; not implemented in this package.

## 学者、软件与本机资料

- Werner Krauth, [Statistical Mechanics: Algorithms and Computations course and examples](https://www.lps.ens.fr/~krauth/index.php/MOOC_SMAC). 本机还有 `Desktop/GitHub/真·Krauth/` 配套示例。小系统精确枚举与逐步验证的组织方式值得学习。
- Synge Todo, [ALPS/looper](https://github.com/wistaria/alps-looper). 这是量子模拟代码，列为软件工程参考，不作为本项目经典算法的实现来源。
- 张潘，[GitHub 主页](https://github.com/panzhang83)。关注其统计推断、张量网络及相关研究；不将其项目笼统视为本包算法实现。
- 邓友金等，[随机团簇研究](https://doi.org/10.1103/PhysRevE.80.036707)。本机 `Desktop/邓友金_dyj/` 有 Worm 讲义与 `CQMC_Ising_program/` 参考程序。前者供算法背景学习，后者提供测量与程序组织的参考；本包实现经过独立验证，并未照搬。
- 本地工作区 `percolation_fss/` 是已有的 Python 项目，包含 Newman–Ziff、binomial convolution 与 FSS 示例，可作为日后交叉验证材料；它不属于本包依赖。
