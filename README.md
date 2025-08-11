# Benchmarking Spike Source Localization Algorithms in High-Density Probes

This repository contains the analysis scripts and supporting functions used in the study:

> Hao Zhao, Xinhe Zhang, Arnau Marin-Llobet, Xinyi Lin, Jia Liu.  
> **Benchmarking Spike Source Localization Algorithms in High-Density Probes**. *PLOS Computational Biology*, 2025.  

## Overview
The code implements and benchmarks three spike source localization algorithms:
- **Center of Mass (COM)**
- **Monopolar Triangulation (MT)**
- **Grid Convolution (GC)**

We evaluate these algorithms using:
1. **Simulated dataset** – generated with [MEArec](https://github.com/alejoe91/MEArec)  
2. **Experimental dataset** – [SPE-1 paired patch-clamp and Neuropixels recordings](https://doi.org/10.1101/370080)

Performance metrics include accuracy, robustness to electrode degradation, and runtime.

## Data
- **Simulated data**: Generated with MEArec (open source) using the included configuration files.  
- **Experimental data**: Download SPE-1 dataset from the [original repository](https://doi.org/10.1101/370080).  
- All required preprocessing scripts are included.

## Usage
Clone this repository and run:
```bash
python run_benchmark.py --config configs/benchmark_config.yaml
