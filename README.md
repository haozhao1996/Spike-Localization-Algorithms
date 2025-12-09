# Spike Localization Algorithms

This repository contains the full pipeline used for our manuscript "Benchmarking spike source localization algorithms in high density probes". We include figure-generation notebooks for both the simulated dataset (MEArec) and experimental dataset (SPE1), as well as a simplified MEArec pipeline with all intermediate recording, spike train, waveform, and localization data directly provided.

---

## Repository Structure

```
Spike-Localization-Algorithms/
│
├── SSL_mearec (2025.12.07).ipynb
├── SSL_mearec (2025.12.07)_AS.ipynb
├── SSL_mearec_toy.ipynb
├── SSL_mearec_toy_AS.ipynb
├── SSL_spe1 (2025.12.07).ipynb
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
└── output/
```

---

# 1. Simulated Dataset Pipeline (MEArec)

Two notebooks are required to run the MEArec benchmarking pipeline:

- **Run pipeline initialized to one set of seeded neuron locations:** `SSL_mearec (2025.12.07).ipynb`
  - Adjust the 'mearec_seed' parameter to an int which seeds the ground-truth neuron locations used in the simulation. In our manuscript, we run this notebook five times with 'mearec_seed' initialized to 42, 43, 44, 45, and 46. 
  - Adjust the 'dead_indices_seed' parameter to a list of ints which seed the pattern of electrode degradation (i.e., the order in which electrodes "drop out"). In our manuscript, we use a list of [42, 43, 44, 45, 46]; each run of the notebook will run trials using five different patterns of electrode degradation.
  - Adjust the 'run_id' parameter to determine the output folder of the run. Do not remove the 'mearec_seed' value appended to the end, which the following file will use to identify localization results.
  - This file generates the simulated dataset, performs all preprocessing and waveform extraction, and performs template and spike localization. Intermediate files (e.g., recording files) are deleted to manage folder size, but localization results are saved.  
- **Aggregating results across multiple trials:** `SSL_mearec (2025.12.07)_AS.ipynb`
  - Once prior pipeline is run with relevant 'mearec_seed' values, populate those values into list in 'mearec_seeds' parameter. In our manuscript, we run this notebook with 'mearec_seeds' set to [42, 43, 44, 45, 46].
  - This files aggregates the localization results across 'mearec_seed' values to generate the benchmarking results presented in the manuscript.

In the interest of presenting an "end to end" dataset including intermediate recording, spike train, and waveform data (which we delete to manage folder size), we also present a simplified toy pipeline:
- **Run pipeline initialized to one set of seeded neuron locations:** `SSL_mearec_toy.ipynb`
  - Identical to `SSL_mearec (2025.12.07).ipynb`, except truncating the degradation scenarios and simulation length.
- **Aggregating results across multiple trials:** `SSL_mearec_toy_AS.ipynb`
  - Identical to `SSL_mearec (2025.12.07)_AS.ipynb`, except truncating the 'mearec_seed' scenarios and simulation length.
- **Intermediate data:** `SSL_mearec_toy_AS.ipynb`
  - We provide an intermediate dataset of recording, spike train, and waveform data generated using simplified pipeline and 'mearec_seed' values of 42 and 43.
  - https://drive.google.com/drive/folders/14EXJ_RqWtAA6alYxuq8WDHz6J-U7QLDh?usp=sharing


# 2. Experimental Dataset Pipeline (SPE-1)

One notebook is required to run the SPE-1 benchmarking pipeline:

- `SSL_spe1 (2025.12.07).ipynb`
  - Data must be uploaded into `spe1/data/` folder. Data for each cell is included in separate folder. We recommend downloading using helpful crcnsget tool (https://github.com/neuromusic/crcnsget).  
  - This file performs all degradation simulation, preprocessing, waveform extraction, and template and spike localization. Intermediate files (e.g., recording files) are deleted to manage folder size, but localization results are saved.  

# 3. Citation

If you use this repository or the simplified MEArec dataset, please cite:

```
Zhao, H., Zhang, X., et al. (2025). Benchmarking spike source localization algorithms in high density probes.
```

---
