# compression_driven_adaptive_sampling
Entropy-driven adaptive sampling for wrist PPG/ACC using an adapted WFPV (Temko, 2017) pipeline (see acknowledgement).

## Contents
- `BlandAltman.m` and the three `PPG_*` scripts come from the WFPV work by Temko, 2017 (used for HR estimation baselines and plots).
- `adaptive_entropy_sampling_last.m` wraps the adapted WFPV code with my finite-state machine (FSM) controller for switching sampling rates based on entropy.
- Entropy helpers: `entropy_proxy_hist.m`, `entropy_proxy_arith_o1.m`, `entropy_proxy_lz.m`.
- Quick sanity checks: `entropy_proxy_sanity_checks.m` (tiny deterministic signals with `assert` checks).

## Data
- IEEE Signal Processing Cup 2015 PPG/ACC dataset (Zhang et al., TROIKA): download from the official sources, e.g., the TROIKA repo and Zenodo mirror:
  - https://github.com/hpi-dhc/TROIKA
  - https://zenodo.org/records/3902710
- Place the MAT files under `Data/` so the `load(['Data/' IDData{...}])` calls work.

## Running the adaptive sampler
1. Ensure the data are in `Data/`.
2. Open `adaptive_entropy_sampling_last.m` and run it in MATLAB. It logs results to `adaptive_entropy_logs.mat` and can generate Bland–Altman plots.

### Switching entropy proxies
- The entropy proxy used for the controller is set at line 119 (`Hacc(i) = entropy_proxy_arith_o1(...)`). Replace with `entropy_proxy_hist` or `entropy_proxy_lz` if desired.

### Tuning adaptive thresholds
- Controller thresholds/hysteresis live at lines 21–28 (`nbits_entropy`, `Th_hi`, `Th_low`, `hi_hold`, `N_look_back`, `N_unstable`, `N_stable`). Adjust to change when the FSM switches between low/high sampling rates.

### Sanity checks
- Run `entropy_proxy_sanity_checks` to print entropy values for simple constant/step/periodic/noisy signals and assert basic ordering/invariance.

## Notes and acknowledgment
- This codebase is heavily commented and cleaned with help from ChatGPT (Codex). Debugging and final design choices are mine; AI assistance was used for refactoring and documentation.
- Baseline WFPV implementation and reference material: https://github.com/andtem2000/PPG (Temko et al., IEEE TBME 2017, “Accurate Heart Rate Monitoring During Physical Exercises Using PPG”). 

## Poster
<img width="823" height="732" alt="image" src="https://github.com/user-attachments/assets/50ce6113-eade-4c91-8389-51dfb96549bb" />
