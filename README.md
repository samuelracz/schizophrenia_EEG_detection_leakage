# Code supplement for "Information Leakage and Performance Overestimation in EEG-Based Schizophrenia Detection: Evidence from Literature and Empirical Analyses"
### Authors: Frigyes Samuel Racz & Gabor Csukly (2025)

This repository contains example codes and corresponding functions for the classification pipelines included in the study "Information Leakage and Performance Overestimation in EEG-Based Schizophrenia Detection: Evidence from Literature and Empirical Analyses". The analyzed EEG datasets are available online and therefore not shared here:
- MSU dataset: Borisov, S. V., et al. ["Analysis of EEG structural synchrony in adolescents with schizophrenic disorders."](https://link.springer.com/article/10.1007/s10747-005-0042-z) Human Physiology 31.3 (2005): 255-261. Data available at: http://brain.bio.msu.ru/eeg_schizophrenia.htm.
- RepOD dataset: Olejarczyk, Elzbieta, and Wojciech Jernajczyk. ["Graph-based analysis of brain connectivity in schizophrenia."](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0188629) PloS one 12.11 (2017): e0188629. Data available at: https://doi.org/10.18150/repod.0107441.
- SU-SZ dataset: Racz, Frigyes Samuel, et al. ["Reduced temporal variability of cortical excitation/inhibition ratio in schizophrenia."](https://www.nature.com/articles/s41537-025-00568-3) Schizophrenia 11.1 (2025): 20. Data available at: https://doi.org/10.5281/zenodo.14808295.

For all three datasets, the following Matlab scripts are provided:
1. Pre-processing of raw EEG data
2. Feature extraction for CNN- and ML-based classification pipelines
3. Scripts performing classification in leaky- and leakage-free implementations

In addition, Matlab workspaces storing ground truth, predicted label and probability score values from the performed evaluations are also provided, which can be used for the reproduction of all results reported in the study.

**Notes:**
- Most scripts require the EEGLAB Toolbox ([Delorme & Makeig, 2004](https://www.sciencedirect.com/science/article/abs/pii/S0165027003003479?casa_token=THNwh9oC20AAAAAA:U0xOJ5c9TUF6Aws36j2sdDQ9f9n3Pv6K8EEgFwLSSlt2RipIDDGY1gGObrahpyGgMtSRjmLxHVeq)) with necessary plug-ins to be added to the matlab path.
- The raw text files of the MSU dataset were already converted to matrices of size 7680*16, and stored as matlab workspaces (expected input for pre-processing script).
