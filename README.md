# Causal brain activity modeling for the understanding of epilepsy pathogenesis

This repository contains the implementation of the causal brain activity modelling described in my bachelor thesis. 

## Installation
```
conda env create -f environment.yml
```

## Requirements for the dataset
1) The analyzed file should be in the European Data Format (EDF)
2) The markings of the IEDs should be recorded in Time Measurement Units (TMu) and save to an event (evt) file. The first column of the file is the corresponding marking.

## Analysis of causal interactions for a single patient 

Sample analysis
```
cd scripts
python run.py --window_length=500 --edf_file=/data/garkots/third_package/pat04_ON_sleep/EEG_1278-export.edf --annotations=/data/garkots/third_package/pat04_ON_sleep/ --path_to_visual_outputs=/home/garkots/invasive_eeg/analysis/visuals/experiment_window_500/pat4/sleep --path_to_results=/home/garkots/invasive_eeg/analysis/results/experiment_window_500/pat4/sleep --path_to_ica_matrix=/home/garkots/invasive_eeg/analysis/matrices_ica/pat4/sleep --state=sleep
```
The following parameters are recommended to be changed from the default values

```
--window_length = 100   # -> length of window for analysis
--path_to_ica_matrix = ...  # -> if an experiment was previously conducted on a patient in order to save computing time it is recommended to specify path to storage of ICA decomposing matrix
--precomputed_ica = True   # -> reads previously computed ICA matrix 
--chosen_components = ICA001 ICA003 ICA023 #-> in case the chosen component are prefered to be chosen manually 
```


Detailed analysis [here](https://docs.google.com/presentation/d/1pNNhWa_rgkRxldGc8owuecJ4eIycxf3pHaG_A0dXX5g/edit?usp=sharing).
