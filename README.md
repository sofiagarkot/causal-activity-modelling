# Causal brain activity modeling for the understanding of epilepsy pathogenesis

This repository contains the implementation of the causal brain activity modelling described in my bachelor thesis. 

## Installation
1. Clone the repository.
2. Create the environment.
```
conda env create -f environment.yml
```
3. Run the analysis with the specified parameters :) 
```
cd scripts
python run.py --window_length=*** --edf_file=***path to edf file*** 
--annotations=***path to directory containing annotations*** 
--path_to_visual_outputs=***path to directory where the results of visuals analysis will be stored (you can turn off the saving of visuals)*** 
--path_to_results=***path to directory where the results (per tau and per component distributions) will be stored*** 
--path_to_ica_matrix=***path to directory where the decomposed components will be stored***
--state=***sleep/wake***
--save_visuals=***True/False***
--automated_choice=***True/False - whther to automatically detect the significant components***
--precomputed_ica=***True/False - in case when the experiment has previously been done on the same patient, to save the time the precomputed decomposed components can be used***
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
