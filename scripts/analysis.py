from patient import Patient
from utils import deviative_choice, per_tau_events_density_plot, get_density_plot, save_to_json

#TODO: implement!
from utils import parse_ica_channels

def make_patient_analysis(args):

    deviative_component_choice = False
    # implement it more carefully
    patient = Patient( edf_file = args.edf_file,
                       annotations = args.annotations,
                       path_to_visual_outputs = args.path_to_visual_outputs,
                       path_to_ica_matrix = args.path_to_ica_matrix,
                       path_to_results = args.path_to_results,
                       window_length = args.window_length,
                       precomputed_ica = args.precomputed_ica,
                       ignore_channels = args.ignore_channels,
                       chosen_components = args.chosen_components,
                       control_experiment = args.control_experiment)


    patient.quantify_annotations()

    patient.run_ica()
    if args.ignore_channels:
        patient.no_channels_during_seizures_states()

    # else:
    #     patient.during_seizures_states()

    patient.normalize_components()


    patient.plot_all_components_during_seizures(savefigures=True)
    if  patient.chosen_components == []:
        if deviative_component_choice:
            chosen = deviative_choice(patient)
        else:
            chosen = parse_ica_channels()
    else:
        chosen = patient.chosen_components

    print("Chosen components names", chosen)

    if args.path_to_visual_outputs:
        patient.plot_selected_components_during_seizures(potentially_significant=chosen, savefigures=True)

    patient.get_significant_components_state(potentially_significant=chosen)
    patient.run_causal_inference()

    tau_max = per_tau_events_density_plot(patient)

    all_density_patient = get_density_plot(patient)

    path_all = patient.path_to_results+'/'+args.state+'_total__.json'

    max_tau_density_patient = get_density_plot(patient, selected_tau=tau_max)

    path_max = patient.path_to_results+'/'+args.state+'_most_frequent_tau__.json'


    save_to_json(path_all, all_density_patient)
    save_to_json(path_max, max_tau_density_patient)
