# Shadow Pricing in Transient Stability-Constrained Optimal Power Flow

Code, input data and results of the paper of the same title
([arXiv:2610.18563](https://arxiv.org/abs/2610.18563)). The repository holds
only what the paper uses.

## Layout

| Folder | Content |
|---|---|
| `TSCOPF/` | The Julia package that builds and solves the models, with the `Manifest.toml` of the environment. |
| `TSCOPF/docs/` | The documentation of the package: the model, equation by equation, and every setting. Published at <https://alexupm95.github.io/tscopf-economic-perspective/dev/>. |
| `kkt/` | Price decomposition: the exports attached to a solve (`kkt_pg_terms.jl`, `kkt_sensitivity.jl`) and the notebooks that turn them into shadow prices and settlements. |
| `INPUT_FILES/` | Input data of the three cases of the paper. |
| `runs/` | The scripts that solve the runs. |
| `results/` | One folder per run printed in the paper. |
| `postprocess/` | The scripts that turn the run folders into the tables and figures of the paper, with their output. |

## Cases and runs

| Case folder | Section | Run | Command |
|---|---|---|---|
| `9bus_linear_H24_6_3` | IV-A | `results/9bus/Results - 2026-10-08 112948` | `julia --project=TSCOPF runs/run_9bus.jl linear` |
| `9bus_colocated_bus1_H24_6_3` | IV-A | `results/9bus/Results - 2026-10-08 113323` | `julia --project=TSCOPF runs/run_9bus.jl colocated` |
| `39bus_gfm` | IV-B | `results/39bus/Results - 2026-10-07 124340` | `julia --project=TSCOPF runs/run_39bus.jl` |

Run the commands from the root of the repository, after
`julia --project=TSCOPF -e "using Pkg; Pkg.instantiate()"` (Julia 1.11). A solve writes a new
`Results - <timestamp>` folder under `new_runs/`, which is not tracked.

Every run folder holds its `run_manifest.toml` (all settings), the input data as solved
(`Inputs/`), the dispatch and its multipliers (`Dispatch/`, with the ACOPF benchmark in
`Dispatch_WarmStart/`), and the trajectories and the multipliers of the dynamic constraints
(`Transient_Stability/`).

`Transient_Stability/dynamic_model_details.txt`, a text dump of the model above 100 MB on
the 39-bus run, is left out for size. A solve rewrites it.

## Tables and figures

From `postprocess/`, with Python 3.11 or later, `pandas` and `matplotlib`:

| Script | Output |
|---|---|
| `tables_9bus.py` | `tables/tab_prices_9bus.tex`, `tables/tab_colocated_bus1_9bus.tex`, `tables/prose_numbers_9bus.txt` |
| `plot_angles_9bus.py` | `figures/angles_tscopf_9bus.pdf` |
| `plot_lambda_9bus.py` | `figures/lambda_colocated_bus1_9bus.pdf` (run `tables_9bus.py` first) |
| `tables_figure_39bus.py` | `tables/tab_prices_39bus_gfm.tex`, `tables/tab_settlement_39bus_gfm.tex`, `figures/angles_psse_39bus_gfm.pdf` |

`psse/aligned.csv` holds the PSS/E time-domain simulation of the 39-bus dispatch, on the time
grid of the TSCOPF run. The figures use the Courier New font.

## Citation

If you use this code or data, please cite the paper:

```bibtex
@misc{coelho2026shadowpricingtransientstabilityconstrained,
      title={Shadow Pricing in Transient Stability-Constrained Optimal Power Flow},
      author={Alex Coelho and Araceli Hernandez and Luis Badesa},
      year={2026},
      eprint={2610.18563},
      archivePrefix={arXiv},
      primaryClass={eess.SY},
      url={https://arxiv.org/abs/2610.18563},
}
```

## Acknowledgements

This work was supported by **MICIU/AEI/10.13039/501100011033** and **ERDF/EU** under grants **PID2023-150401OA-C22** and **PID2022-141609OB-I00**, as well as by the Madrid Government (Comunidad de Madrid-Spain) under the Multiannual Agreement 2023-2026 with Universidad Politécnica de Madrid, `Line A - Emerging PIs' (grant number: **24-DWGG5L-33-SMHGZ1**).

<div align="center" style="display:flex; flex-wrap:wrap; align-items:center; justify-content:center; gap:24px;">
  <img src="./utils/logo_MICIU_AEI.png" alt="MICIU/AEI" height="200" style="object-fit:contain;">
  <img src="./utils/logo_CM.png" alt="Comunidad de Madrid" height="200" style="object-fit:contain;">
</div>
