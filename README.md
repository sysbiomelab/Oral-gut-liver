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
│   ├── Figure_1.Rmd            R markdown for For Figure 1 visualisations
│   ├── dmngut.RData             # Gut DMN model
│   ├── dmnoral.RData            # Oral DMN model
│   ├── GEM_modelling.m          # Single-species GEM simulations
│   └── Tables/                  # Data for Figure 1
│
├── Figure_2/
│   ├── Figure_2.Rmd           R markdown fo# For Figure 2 visualisations
│
│   ├── GSMM/
│   │   ├── Bacterial_community_generate.m   # Community model construction
│   │   ├── run_host.m                       # Host GSMM simulations
│   │   ├── data.mat                         # Model input data
│   │   ├── RPMI_humanGEM.csv                # Liver GSMM constrains
│   │   ├── eflux_normal.csv                 # Flux (healthy)
│   │   ├── eflux_cirrhosis.csv              # Flux (cirrhosis)
│   │   └── Average_Low_SampleID.xlsx        # Input data for community GSMM construction
│   │
│   └── Tables/              # Data for Figure 2
```

---

## Requirements

### R

* `DirichletMultinomial`
* `ComplexHeatmap`
* `ggplot2`

### MATLAB

* MATLAB (R2023 or later recommended)
* COBRA Toolbox
* MIGRENE Toolbox

---

## Data Availability

All processed data required to reproduce the figures are included in:

* `Figure_1/Tables/`
* `Figure_2/Tables/`

---

## Citation

If you use this code, please cite:

> *Citation will be added upon publication*

---

## Contact

**Yi Jin**
King’s College London  
PhD student in Computational and Systems Biology

---

## License

This repository is provided for academic research use. For other uses, please contact the authors.

---

