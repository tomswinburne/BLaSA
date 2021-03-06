## BLaSA - Bond Lattice Sampling and Analytical Approximation
:copyright: TD Swinburne 2020 MIT License swinburne@cinam.univ-mrs.fr

Execution and analysis of bond lattice dynamics. When using this code, please cite

> *Anharmonic free energy of lattice vibrations in fcc crystals from a mean field bond*   
> Accepted for publication in PRB Rapid Communications, August 2020
> Thomas D Swinburne, Jan Janssen, Mira Todorova, Gideon Simpson, Petr Plechac, Mitchell Luskin and Jörg Neugebauer


## Setting the parameters
- All parameters are in `input.json`. Default parameters are those used in the above publication.
- Energy units are electron volts, length units are Angstrom.

- `dump_folder` : output data will be stored in `run/dump_folder`, which will be created if required
- `workers_per_value` : number of workers assigned to each parameter value. Must factorize total number of cores
- `md` -> `bins` :  number of histogram bins
- `md` -> `n` : supercell has `n^3` unit cells
- `md` -> `therm`, `steps`   : number of timesteps to thermalize and sample
- `potential` -> `D0`, `a0`, `AL` : Morse potential `V(r) = D0 * exp(-AL*(r-a0))`
- `potential` -> `r_min`, `r_max` : limits for histogram
- `am_array`: list of strains to sample. `am=0.01` is 1% tensile strain
- `t_array`: list of temperatures
- `RT_array` : list of transverse strength ratios to sample.
- `JointHist` : 0/1. Turns on recording of joint histograms to investigate bond-bond correlations. Returns three `bins*bins` arrays.

## Cloud hosted analysis notebooks
- This repository is hosted online using binder - [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tomswinburne/BLaSA/master?filepath=analyze)
- Due to the CPU hour restrictions of this service we do not recommend generating high quality reference data online, only running the analysis notebooks

## Automated data generation and analysis with Snakemake (recommended)
- Requires a working [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) package manager
- Install Snakemake workflow manager
```
conda install -c bioconda -c conda-forge snakemake
```
- Run calculations with
```
snakemake --use-conda --cores 10
```
where `--cores` is the total number of CPUs. This must be an integer multiple of `workers_per_value`

## Manual data generation and analysis
- Instructions can be found [here](ManualInstallation.md)
