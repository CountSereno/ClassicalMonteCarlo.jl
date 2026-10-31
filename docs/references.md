# 参考资料与来源

实现根据公开的算法定义独立编写；参考资料不意味着作者认可本项目。论文、示例代码分别受其各自许可证及版权条款约束。

## 经典算法与验证

- Metropolis et al., “Equation of State Calculations by Fast Computing Machines” (1953), [DOI: 10.1063/1.1699114](https://doi.org/10.1063/1.1699114).
- Swendsen and Wang, “Nonuniversal critical dynamics in Monte Carlo simulations” (1987), [DOI: 10.1103/PhysRevLett.58.86](https://doi.org/10.1103/PhysRevLett.58.86).
- Wolff, “Collective Monte Carlo Updating for Spin Systems” (1989), [DOI: 10.1103/PhysRevLett.62.361](https://doi.org/10.1103/PhysRevLett.62.361).
- Sweeny, “Monte Carlo study of weighted percolation clusters relevant to the Potts models” (1983), [DOI: 10.1103/PhysRevB.27.4445](https://doi.org/10.1103/PhysRevB.27.4445).
- Beffara and Duminil-Copin, “The self-dual point of the two-dimensional random-cluster model is critical for q≥1”, [arXiv:1006.5073](https://arxiv.org/abs/1006.5073). 多模型实验代入方格的已知临界参数，未重新拟合临界点。
- Elçi and Weigel, “Efficient simulation of the random-cluster model” (2013), [arXiv:1307.6647](https://arxiv.org/abs/1307.6647). 参考 FK 表示、几何观测量和连通性成本的区分；本包独立实现双端交替 BFS，但尚未实现该工作的动态森林算法。
- Xu, Salas and Deng, “Finite-Size Corrections from the Subleading Magnetic Scaling Field for the Ising and Potts Models in Two Dimensions” (2025), [DOI: 10.3390/e27040418](https://doi.org/10.3390/e27040418). 团簇分形维数及多种修正项背景；该论文的 loop 模型修正振幅不直接用于本项目方格随机团簇模型。
- Deng, Garoni and Sokal, “Critical speeding-up in a local dynamics for the random-cluster model” (2007), [arXiv:cond-mat/0701113](https://arxiv.org/abs/cond-mat/0701113). 本轮每 sweep 测量一次，不能分辨其研究中的子 sweep 快过程，亦未测定动力学指数。
- Prokof’ev and Svistunov, “Worm Algorithms for Classical Statistical Models” (2001), [arXiv](https://arxiv.org/abs/cond-mat/0103146).
- Wang and Landau, “Determining the density of states for classical statistical models” (2001), [arXiv](https://arxiv.org/abs/cond-mat/0107006).
- Newman and Ziff, [percolation algorithm and original implementation](https://websites.umich.edu/~mejn/percolation/). A useful future implementation and validation target; not implemented in this package.

