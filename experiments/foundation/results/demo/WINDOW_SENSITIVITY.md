# Offline ACF-window sensitivity at L=24

Reanalysed 12 saved critical-temperature chains with maximum lags 512, 1024, 2048. No samples were generated and the primary analysis remains at lag 512. All baseline tau values were reproduced exactly; input and analysis-source hashes are recorded in analysis/window_sensitivity.toml.

| Algorithm | Observable | Mean tau at 512 | At 1024 | At 2048 | Truncated at 2048 |
|---|---|---:|---:|---:|---:|
| metropolis | energy | 49.422 | 50.152 | 50.152 | 0 |
| metropolis | absM | 115.681 | 117.957 | 117.957 | 0 |
| wolff | energy | 3.334 | 3.334 | 3.334 | 0 |
| wolff | absM | 2.678 | 2.678 | 2.678 | 0 |
| swendsen_wang | energy | 3.631 | 3.631 | 3.631 | 0 |
| swendsen_wang | absM | 3.293 | 3.293 | 3.293 | 0 |

Tau uses native algorithm updates, so different algorithms have different clocks. For a fixed chain the inferred ESS varies inversely with tau; timing data are reused without retiming.

| Changed chain | Observable | Tau ratio 2048 / 512 | ESS ratio | Last included lag |
|---|---|---:|---:|---:|
| run_038 (metropolis) | energy | 1.0487 | 0.9535 | 679 |
| run_038 (metropolis) | absM | 1.0629 | 0.9408 | 794 |
| run_039 (metropolis) | energy | 1.0008 | 0.9992 | 520 |

This follow-up tests one window choice on the same finite traces. Reaching a first zero before the new limit does not establish stationarity, remove long-tail noise, or validate signed-M susceptibility errors. It does not replace longer independent chains or a separate equilibration study.
