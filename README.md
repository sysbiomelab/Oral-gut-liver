# Oral–Gut–Liver Microbiome Modelling

Code accompanying the study:

**“Integrative host–microbiome modelling uncovers the implication of oral–gut translocation in advanced cirrhosis”**

---

## Overview

This repository provides the analysis workflows and modelling code used to investigate oral–gut microbiome translocation and estimate its metabolic implications in liver cirrhosis.

The study integrates multi-scale analyses, including:

* Reactobiome profiling and reactotyping
* Oral–gut metabolic distance (OGMD) analysis
* Genome-scale metabolic modelling (GEMs) of microbial species
* Community-level metabolic simulations
* Host–microbiome metabolic interaction simulations (liver, brain, muscle)

---

## Repository Structure

```bash
Oral-gut-liver/
├── Figure_1/
│   ├── Figure_1.Rmd             # R markdown for Figure 1 visualisations
│   ├── dmngut.RData             # Gut DMN model
│   ├── dmnoral.RData            # Oral DMN model
│   ├── GEM_modelling.m          # Single-species GEM simulations
│   └── Tables/                  # Data for Figure 1
│
├── Figure_2/
│   ├── Figure_2.Rmd             # R markdown for Figure 2 visualisations
│
│   ├── GSMM/
│   │   ├── Bacterial_community_generate.m   # Community model construction
│   │   ├── run_host.m                       # Host GSMM simulations
│   │   ├── data.mat                         # Model input data
│   │   ├── RPMI_humanGEM.csv                # Liver GSMM constraints
│   │   ├── eflux_normal.csv                 # Flux constraints (healthy)
│   │   ├── eflux_cirrhosis.csv              # Flux constraints (cirrhosis)
│   │   └── Average_Low_SampleID.xlsx        # Input data for community GSMM construction
│   │
│   └── Tables/                              # Data for Figure 2
```

---

## Requirements

### R

Required R packages:

* `DirichletMultinomial`
* `ComplexHeatmap`
* `ggplot2`
* `vegan`
* `tidyverse`

### MATLAB

* MATLAB (R2023a or later recommended)
* COBRA Toolbox
* MIGRENE Toolbox

---

## Reproducibility

All analyses were performed using custom R and MATLAB scripts as described in the manuscript and supplementary methods.

All processed data required to reproduce the figures are included in:

* `Figure_1/Tables/`
* `Figure_2/Tables/`

The repository contains scripts for:

* Reactobiome construction and reactotyping
* Oral–gut metabolic distance (OGMD) analysis
* Single-species GEM simulations
* Community-level GSMM construction
* Host metabolic simulations under ammonia and acetate perturbations

---

## Data Availability

The shotgun metagenomic datasets analysed in this study are publicly available from the European Nucleotide Archive (ENA):

* GLA cohort: PRJEB52891
* RIFSYS cohort: PRJEB38481
* Healthy cohort: PRJEB38483

The corresponding genome-scale metabolic models (GSMMs) are available from:

https://www.microbiomeatlas.org

---

## Citation

If you use this repository in your work, please cite:

> Jin Y, et al. 2026 *Integrative host–microbiome modelling uncovers the implication of oral–gut translocation in advanced cirrhosis*. iMeta e70131. https://doi.org/10.1002/imt2.70131

Online article:

https://onlinelibrary.wiley.com/doi/full/10.1002/imt2.70131

---

## Contact

**Yi Jin**  
Centre for Host–Microbiome Interactions  
Faculty of Dentistry, Oral & Craniofacial Sciences  
King’s College London  
PhD student in Computational and Systems Biology

---

## License

This repository is provided for academic research use only.

For commercial use or collaboration enquiries, please contact the authors.

---
