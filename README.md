# Assessing Careless and Insufficient Effort Responding

This repository contains the supplementary files for the manuscript: 
> Identifying Careless Responding in Ecological Momentary Assessment: Inconsistent Signals from Different Detection Methods in the WARN-D Data


You can cite this project as: 

```BibTeX
Will be updated once preprint is out. 
```

## File structure


### `\data\`
Contains the data used for this project. As described in the manuscript, data used for this project wil be made available once all data collection is finished and all data are cleaned and anonymized. 


### `\scripts\`

All Quarto files are available as source code (`.qmd`) and as a rendered version (`.html`)

| Script                          | Description                                                                 |
|---------------------------------|-----------------------------------------------------------------------------|
| `00_functions.R`                | Contains auxiliary functions                                                |
| `01_calculate_indices`          | Code to calculate indices of C\\IER                                         |
| `01_visualize_indices`          | Code to analyze and visualize all indices and model results. Produces main figures for the manuscript |
| `02_LPA`                        | Code and visualizations for the latent profile analysis                     |
| `03_Mixture_IRT_model`          | Code and visualizations for the mixture IRT model                           |
| `04_Response_Time_Mixture_Model`| Code and visualizations for the response time mixture model                 |
| `screentimeCIER.stan`           | Stan code for the response time mixture model                               |


### `\figures\`

Contains all figures used in the manuscript and in the supplementary code. 