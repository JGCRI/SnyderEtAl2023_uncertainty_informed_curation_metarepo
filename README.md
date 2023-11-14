[![DOI](https://zenodo.org/badge/478588973.svg)](https://zenodo.org/badge/latestdoi/478588973)

# SnyderEtAl2023_uncertainty_informed_curation_metarepo
 
**Uncertainty-informed selection of CMIP6 Earth System Model subsets for use in multisectoral and impact models**

Abigail Snyder<sup>1\*</sup>,  Noah Prime<sup>1</sup>, Claudia Tebaldi<sup>1</sup>, and Kalyn Dorheim<sup>1</sup>

<sup>1 </sup>  Joint Global Change Research Institute, Pacific Northwest National Laboratory and University of Maryland, College Park, MD

\* corresponding author:  abigail.snyder@pnnl.gov


## Abstract
Earth System Models (ESMs)  are heavily used to provide inputs to specific impact and multisectoral dynamic models. Therefore, representing the full range of model uncertainty, scenario uncertainty, and interannual variability that ensembles of  ESMs capture ,  is critical to the exploration of the future co-evolution of the integrated human-Earth system. The pre-eminent source of these ensembles has been the Coupled Model Intercomparison Project (CMIP) however, the increased participation by modeling centers in every new phase of CMIP relative to previous eras  results in increasingly large archives of ESM data that can be intractable for impact modelers to effectively utilize, due to computational constraints and the challenges of analyzing large datasets. In this work, we describe a method to select a subset of the latest phase, CMIP6, models for use as inputs to a sectoral impact or multisectoral model, while still  representing the range of  model uncertainty, scenario uncertainty, and interannual variability of the full CMIP6 ESM results. This method is intended to help human-relevant impact and multisectoral modelers select climate information from the CMIP archive efficiently. This is particularly critical for large ensemble experiments of multisectoral dynamic models that may be varying additional features beyond climate inputs in a factorial design, thus putting constraints on the number of climate simulations that can be used. We focus on temperature and precipitation outputs of ESMs, as these are two of the most used variables among impact models and many other key input variables for impacts are at least correlated with one or both of temperature and precipitation (e.g. relative humidity). Besides preserving the multi-model ensemble variance characteristics, we prioritize selecting ESMs in the subset that preserve the very likely distribution of equilibrium climate sensitivity values as assessed by the latest IPCC report. This approach could be applied to other output variables of ESMs and, when combined with emulators, offers a flexible framework for designing more efficient experiments on human-relevant climate impacts. It can also provide greater insight into the properties of existing ESMs and the method may be informative for future experiment planning across ESMs. 

This repository collects code for the calculations performed in this paper.

## Journal reference
Submitted to Earth System Dynamics 


## Workflow

Run the following scripts to re-create the manuscript experiments:

| Script Name | Description | 
| --- | --- | 
| `aoi_extractions.py` | Script that processes NetCDF gridded time series data of multiple variables(T, P) to time series of annual average values in each of the IPCC WGI regions. | 
| `metric_calc.ipynb` | Script that processes annual time series data in each IPCC region to for  | 
| `A.intermediate_exp.py` | Script that emulates intermediate scenarios using ssp126 and ssp585 runs as the archive, GSAT outputs  | 
| `AB.experiment_tolerance_sweep.py` | Script that emulates intermediate scenarios using all available CMIP6 runs over a range of tolerances and calculates error statistics, GSAT outputs and summary statistics outputs  | 
| `B.ICEnsembles_ESD.r` | Script evaluating the outputs of `A.inital_cond_exp.py`, calculating error statistics and plotting results | 
| `B.IntermediateScenarios_GSAT_ESD.r` | Script evaluating the outputs of `A.intermediate_exp.py`, calculating error statistics and plotting results | 
| `B.IntermediateScenarios_Gridded_ESD.r` | Script evaluating the outputs of `A.tas_psl_pr.py`, calculating SOI and error statistics and plotting results | 


Data generated from the A and AB steps for this publication are archived (https://doi.org/10.5281/zenodo.6461693).


