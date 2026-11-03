# tools/ and results/

Scripts used for the paper experiments (MATLAB/Octave; run from the repository root).

| script | purpose |
|---|---|
| `run_cell.m` | runs one (spacing, deployment) Monte Carlo cell and saves `results/cells/cell_<tag>.mat` (per-epoch detector statistics, mission ids). Octave CLI: `octave run_cell.m <g> <Nruns> <alt> <circle|los> <R_los> <alhourani|freespace> <traj_scale> <Pt_spoof or 0> <tag>` |
| `aggregate_spacing.m` | Pd (+95% bootstrap CI over missions), AUC, distortion for the spacing sweep cells (`cell_g*.mat`); 7 detectors incl. 1-min windowed variants |
| `aggregate_sweep.m` | same for the transmitter-altitude / LOS-placement cells |
| `validate_dilation.m` | first-order dilation model vs measured distortion on a noiseless mission (Fig. "validation") and the noiseless scenario dump |
| `sweep_altitude.m` | noiseless dilation (model + measured) vs transmitter altitude and for the LOS placement |
| `diag_deployment.m` | shows that the distortion of a static 2.5 km array correlates with the swarm's distance from the array centre (inside/outside) |
| `verify_cell_50m.m` | reproduces mission 1 of the 50 m cell with the current code and checks the saved statistics bit for bit |

`results/` holds the aggregated CSVs used for the figures and the Monte Carlo cells (`results/cells/*.mat`, ~0.5 MB each).
