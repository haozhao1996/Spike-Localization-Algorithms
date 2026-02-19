# Spike Localization Algorithms

This repository contains the full pipeline used for our manuscript "Benchmarking spike source localization algorithms in high density probes". We include figure-generation notebooks for both the simulated ground truth dataset (MEArec) and experimental ground truth dataset (SPE1), along with intermediate spike train, spike and template waveform, and localization data.

---

## Repository Structure

```
Spike-Localization-Algorithms/
│
├── SSL_mearec (2026.02.10).ipynb
├── SSL_mearec (2026.02.10)_AS.ipynb
├── SSL_spe1 (2026.02.10).ipynb
│
├── src/
│   └── util_eval.py
├── spe1/
│   ├── data/
│       └── c14/
│       └── c15/
│       └── .../
│   └── Data_Summary.xlsx
│   └── chanMap.mat
│
├── output/
│   ├── 20260210/
│       └── spe1/
│   └── 20260210_42/
│   └── 20260210_43/
│   └── .../
```
---

# 1. Simulated Dataset Pipeline (MEArec)

Two notebooks are required to run the MEArec benchmarking pipeline:

- **Run pipeline initialized to one set of seeded neuron locations:** `SSL_mearec (2026.02.10).ipynb`
  - Adjust the 'mearec_seed' parameter to an int which seeds the ground-truth neuron locations used in the simulation. In our manuscript, we run this notebook five times with 'mearec_seed' initialized to 42, 43, 44, 45, and 46. 
  - Adjust the 'dead_indices_seed' parameter to a list of ints which seed the pattern of electrode degradation (i.e., the order in which electrodes "drop out"). In our manuscript, we use a list of [42, 43, 44, 45, 46]; each run of the notebook will run trials using five different patterns of electrode degradation.
  - Adjust the 'run_id' parameter to determine the output folder of the run. Do not remove the 'mearec_seed' value appended to the end, which the following file will use to aggregate localization results across mearec_seed values.
  - This file generates the simulated dataset, performs all preprocessing and waveform extraction, and performs template and spike localization. The majority of recording files are deleted to manage folder size, but intermediate results including spike trains, waveform data, and localization results are saved.  
- **Aggregating results across multiple trials:** `SSL_mearec (2026.02.10)_AS.ipynb`
  - Once prior pipeline is run with relevant 'mearec_seed' values, populate those values into list in 'mearec_seeds' parameter. In our manuscript, we run this notebook with 'mearec_seeds' set to [42, 43, 44, 45, 46].
  - This files aggregates the localization results across 'mearec_seed' values to generate the benchmarking results presented in the manuscript.

In the interest of presenting an "end to end" dataset, we provide an intermediate dataset of spike train, spike and template waveform, and localization data (https://drive.google.com/drive/folders/14EXJ_RqWtAA6alYxuq8WDHz6J-U7QLDh?usp=sharing). The MEArec data is saved in `20260210_42`, which represents the intermediate data using mearec_seed=42. We would run the pipeline with the four other mearec_seed values to produce four other intermediate folders (e.g., `20260210_43`, `20260210_44`, etc.), and then run `SSL_mearec (2026.02.10)_AS.ipynb` which aggregates localization results and save figures in `20260210`.  

# 2. Experimental Dataset Pipeline (SPE-1)

One notebook is required to run the SPE-1 benchmarking pipeline:

- `SSL_spe1 (2026.02.10).ipynb`
  - Data must be uploaded into `spe1/data/` folder. Data for each cell is included in separate folder. We recommend downloading using helpful crcnsget tool (https://github.com/neuromusic/crcnsget).  
  - This file performs all degradation simulation, preprocessing, waveform extraction, and template and spike localization. The majority of recording files are deleted to manage folder size, but intermediate results including spike trains, waveform data, and localization results are saved.

In the interest of presenting an "end to end" dataset, we provide an intermediate dataset of spike train, spike and template waveform, and localization data (https://drive.google.com/drive/folders/14EXJ_RqWtAA6alYxuq8WDHz6J-U7QLDh?usp=sharing). The SPE-1 data is saved in `20260210/spe1`.

# 3. Citation

If you use this repository or intermediate datasets, please cite:

```
Zhao, H., Zhang, X., et al. (2025). Benchmarking spike source localization algorithms in high density probes.
```

---
