import mne
import numpy as np
import pandas as pd
from tqdm import tqdm
import random

import matplotlib.pyplot as plt

from pingouin import partial_corr, corr
from sympy.utilities.iterables import multiset_permutations

from utils import translate, get_average, no_channels_read_annotations, save_to_json
from errors import OutOfSizeError

class Patient:


    def __init__(self, edf_file,
                 annotations,
                 path_to_visual_outputs,
                 path_to_ica_matrix,
                 path_to_results,
                 window_length,
                 ignore_channels,
                 precomputed_ica,
                 chosen_components,
                 control_experiment):

        self.SEED = 97

        self.chosen_components = chosen_components
        self.control_experiment = control_experiment

        print("CHOSEN COMPONENTS", self.chosen_components)

        self.window_to_detect = window_length
        self.precomputed_ica = precomputed_ica

        self.raw = mne.io.read_raw_edf(edf_file)
        self.raw_data = self.raw.get_data()
        # .filter(l_freq=1., h_freq=None)
        self.fs = self.raw.info['sfreq']
        self.ch_names = self.raw.ch_names

        self.annotations = no_channels_read_annotations(annotations,
                                                            fs = self.fs,
                                                            window=window_length,
                                                            max_length=self.raw_data.shape[1])

        if control_experiment:
            key_random_annotations = [i for i in self.annotations.keys()][0]
            len_values = len(self.annotations[key_random_annotations])
            new_vals = [random.randrange(window_length, self.raw_data.shape[1]-window_length) for i in range(len_values)]
            self.annotations[key_random_annotations] = new_vals

        self.path_to_visual_outputs = path_to_visual_outputs
        self.path_to_results = path_to_results
        self.path_to_ica_matrix = path_to_ica_matrix



    def quantify_annotations(self):

        for k in self.annotations:
            print('Channel ', k, 'with ', len(self.annotations[k]), 'markings')

        self.all_markings = set()
        for timestamps in self.annotations:
            self.all_markings = self.all_markings.union(set(self.annotations[timestamps]))
        print("Total number of markings", len(self.all_markings))


        general_info = {'fs':self.fs,
                        'number of markings':len(self.all_markings)}

        save_to_json(self.path_to_results+'/general_info.json', general_info)

    def run_ica(self, variability_explained=0.95):
        '''
        Function that runs the ICA algorithm on the raw data. Saves the decomposed sources into the variable self.ica_sources_data
        :param variability_explained: float
        '''
        if self.precomputed_ica:
            self.ica_sources_data = np.load(self.path_to_ica_matrix+'/ica_matrix.npy')
            self.topomap = np.load(self.path_to_ica_matrix+'/ica_topomap.npy')

            self.ica_sources_ch_names = []
            with open(self.path_to_ica_matrix+"/component_names.txt", "r") as f:
                for line in f:
                    self.ica_sources_ch_names.append(line.strip())

        else:
            ica = mne.preprocessing.ICA(n_components=variability_explained, random_state=self.SEED)
            ica.fit(self.raw)

            self.topomap = ica.get_components()

            ica_sources = ica.get_sources(self.raw)
            self.ica_sources_ch_names = ica_sources.ch_names
            self.ica_sources_data = ica_sources.get_data()

            np.save(self.path_to_ica_matrix+'/ica_matrix.npy', self.ica_sources_data)
            np.save(self.path_to_ica_matrix+'/ica_topomap.npy', self.topomap)

            with open(self.path_to_ica_matrix+"/component_names.txt", "w") as f:
                for s in self.ica_sources_ch_names:
                    f.write(str(s) + "\n")

            self.precomputed_ica = True

    def no_channels_during_seizures_states(self):
        '''
        The function is used to derive and save the patient's state during IEDs into  a dictionary: self.per_component_state_during_seizures
        Function stores the state of all channels and componets before and after IED in scope of a window
        Created for analysis of channels with no prior text-processing
        :return:
        '''

        per_channel_ieds, per_component_state_during_seizures = {}, {}
        for i, channel in enumerate(self.ch_names):

            data = self.raw_data[i]
            for timestamp in self.all_markings:
                pointer = timestamp

                ranged_data = data[pointer - self.window_to_detect:pointer + self.window_to_detect]

                if len(ranged_data) != self.window_to_detect * 2:
                    raise OutOfSizeError(pointer)

                if channel not in per_channel_ieds.keys():
                    per_channel_ieds[channel] = [ranged_data]
                else:
                    per_channel_ieds[channel].append(ranged_data)

        self.per_channel_ieds = per_channel_ieds
        if self.precomputed_ica:

            for i, ch in enumerate(self.ica_sources_ch_names):

                data = self.ica_sources_data[i]

                for timestamp in self.all_markings:

                    pointer = timestamp

                    ranged_data = data[pointer - self.window_to_detect:pointer + self.window_to_detect]

                    if ch not in per_component_state_during_seizures.keys():
                        per_component_state_during_seizures[ch] = [ranged_data]
                    else:
                        per_component_state_during_seizures[ch].append(ranged_data)

            self.per_channel_ieds, self.per_component_state_during_seizures = per_channel_ieds, per_component_state_during_seizures

    # def during_ieds_states(self):
    #     '''
    #     The function is used to derive and save the patient's state during IEDs into  a dictionary: self.per_component_state_during_seizures
    #     :return:
    #     '''
    #     per_channel_ieds, per_component_state_during_seizures = {}, {}
    #
    #     print("Text-processed channels:")
    #     # Some channels are not text-processed conviniently and need the following text-processing
    #     for channel in self.annotations:
    #         print(channel)
    #
    #         try:
    #             for ch_raw in self.ch_names:
    #                 if channel.upper() in ch_raw and channel not in per_channel_ieds.keys():
    #                     channel = ch_raw
    #
    #             data = self.raw_data[self.ch_names.index(channel.upper())]
    #         except:
    #             continue
    #
    #         for timestamp in self.annotations[channel[:3]]:
    #
    #             pointer = int(timestamp * self.fs)
    #
    #             ranged_data = data[pointer - self.window_to_detect:pointer + self.window_to_detect]
    #
    #             if channel not in per_channel_ieds.keys():
    #                 per_channel_ieds[channel] = [ranged_data]
    #             else:
    #                 per_channel_ieds[channel].append(ranged_data)
    #
    #     for i, ch in enumerate(self.ica_sources_ch_names):
    #
    #         data = self.ica_sources_data[i]
    #
    #         for timestamp in self.all_markings:
    #             pointer = int(timestamp * self.fs)
    #
    #             ranged_data = data[pointer - self.window_to_detect:pointer + self.window_to_detect]
    #
    #             if ch not in per_component_state_during_seizures.keys():
    #                 per_component_state_during_seizures[ch] = [ranged_data]
    #             else:
    #                 per_component_state_during_seizures[ch].append(ranged_data)
    #
    #     self.per_channel_ieds, self.per_component_state_during_seizures = per_channel_ieds, per_component_state_during_seizures

    def normalize_components(self):
        '''
        The function is used to normalize the components, s.t. while visualizing they will be comparable to the signal on the electrodes during IEDs.
        '''
        normalized_per_component_states = {}
        leftMin, leftMax = np.min(self.ica_sources_data), np.max(self.ica_sources_data)
        rightMin, rightMax = np.min(self.raw_data), np.max(self.raw_data)

        for c in self.per_component_state_during_seizures:
            result = []

            for seizure in self.per_component_state_during_seizures[c]:
                normalized_seizure = []
                for el in seizure:
                    normalized_seizure.append(translate(el, leftMin, leftMax, rightMin, rightMax))
                result.append(normalized_seizure)
            normalized_per_component_states[c] = np.array(result)

        self.normalized_per_component_states = normalized_per_component_states

    def plot_all_components_during_seizures(self, savefigures=False):
        '''
        The function is used to plot the states of the components during all the IEDs.
        '''

        for channel in self.per_channel_ieds:

            average_signal, signal_stdev = get_average(self.per_channel_ieds[channel])

            duration_timestamp = [i - self.window_to_detect for i in range(len(average_signal))]

            fig = plt.figure(figsize=(24, 10))

            plt.plot(duration_timestamp, average_signal, lw=3, label='pattern')
            lower_bound = average_signal - signal_stdev
            upper_bound = average_signal + signal_stdev
            plt.fill_between(duration_timestamp, lower_bound, upper_bound, facecolor='lightgoldenrodyellow', alpha=0.5,
                             label='1 sigma range')

            plt.title("Average state during IEDs on channel " + str(channel), fontsize=45)

            for component in self.normalized_per_component_states:
                average_component, component_stdev = get_average(self.normalized_per_component_states[component])

                plt.plot(duration_timestamp, average_component, lw=1, label=component)

            plt.title("Average state during IEDs of channel " + str(channel), fontsize=45)
            plt.legend(prop={'size': 7})
            plt.xlabel('Time relative to an IED', fontsize=40)
            plt.ylabel('Signal, mV', fontsize=40)

            plt.xticks(fontsize = 30)
            plt.yticks(fontsize=30)

            if savefigures:
                plt.savefig(self.path_to_visual_outputs+'/all_components/ch_'+str(channel)+'.jpg')
            plt.close()

    def plot_selected_components_during_seizures(self, potentially_significant, savefigures = False):
        '''
        The function plots channels stated as well as states on selected components during IEDs.
        :param potentially_significant: list(str)
        '''
        for channel in self.per_channel_ieds:

            average_signal, signal_stdev = get_average(self.per_channel_ieds[channel])
            duration_timestamp = [i for i in range(len(average_signal))]

            fig = plt.figure(figsize=(24, 10))

            plt.plot(duration_timestamp, average_signal, lw=3, label='pattern')
            lower_bound = average_signal - signal_stdev
            upper_bound = average_signal + signal_stdev
            plt.fill_between(duration_timestamp, lower_bound, upper_bound, facecolor='lightgoldenrodyellow', alpha=0.5,
                             label='1 sigma range')

            for component in self.normalized_per_component_states:

                if component in potentially_significant:
                    average_component, component_stdev = get_average(self.normalized_per_component_states[component])
                    plt.plot(duration_timestamp, average_component, lw=1, label=component)

            plt.title("Average state during IEDs on channel " + str(channel), fontsize=45)
            plt.legend()
            plt.xlabel('Time relative to an IED', fontsize=40)
            plt.ylabel('Signal, mV', fontsize=40)

            plt.xticks(fontsize = 30)
            plt.yticks(fontsize=30)
            plt.legend(prop={'size': 30})

            if savefigures:
                plt.savefig(self.path_to_visual_outputs+'/significant_components/ch_'+str(channel)+'.jpg')
            plt.close()


    def plot_components_during_seizures_nochannels(self, potentially_significant=False, savefigures = False):
        '''
        The function plots only selected components during IEDs.
        :param potentially_significant: list(str)
        '''
        if potentially_significant:
            fig = plt.figure(figsize=(24, 10))
            plt.title('Average components state during seizures', fontsize=23)
            for component in self.per_component_state_during_seizures:
                if component in potentially_significant:
                    average_component, component_stdev = get_average(
                        self.normalized_per_component_states[component])
                    duration_timestamp = [i for i in range(len(average_component))]
                    plt.plot(duration_timestamp, average_component, lw=1, label=component)

            if savefigures:
                plt.savefig(self.path_to_visual_outputs + '/potentially_significant_during_seizures.jpg')
            plt.close()

        else:
            fig = plt.figure(figsize=(24, 10))
            plt.title('Average components state during seizures', fontsize=23)
            for component in self.per_component_state_during_seizures:
                average_component, component_stdev = get_average(self.normalized_per_component_states[component])
                duration_timestamp = [i for i in range(len(average_component))]
                plt.plot(duration_timestamp, average_component, lw=1, label=component)

            if savefigures:
                plt.savefig(self.path_to_visual_outputs + '/all_components_during_seizures.jpg')
            plt.close()

    def plot_on_channel(self, ch_index_list, savefigures = False):

        fig = plt.figure(figsize=(24, 10))

        for i, channel in enumerate(self.per_channel_ieds):

            if i in ch_index_list or channel in ch_index_list:
                print("Channel ", channel)

                average_signal, signal_stdev = get_average(self.per_channel_ieds[channel])

                duration_timestamp = [i - self.window_to_detect for i in range(len(average_signal))]

                fig = plt.figure(figsize=(24, 10))

                plt.plot(duration_timestamp, average_signal, lw=3, label='pattern')
                lower_bound = average_signal - signal_stdev
                upper_bound = average_signal + signal_stdev
                plt.fill_between(duration_timestamp, lower_bound, upper_bound, facecolor='lightyellow', alpha=0.5,
                                 label='1 sigma range')

                plt.title("Channel " + str(channel), fontsize=23)

        if savefigures:
            plt.savefig(self.path_to_visual_outputs + '/selected_channels'+str(ch_index_list)+'.jpg')
        plt.close()

    def get_significant_components_state(self, potentially_significant, to_save = False):
            '''
            After the significant components were manually chosen,
            this function saves them into one variable stored withing an object of a class.
            :param potentially_significant: list(str)
            '''
            self.significant_components = {}

            for component in self.normalized_per_component_states:
                if component in potentially_significant:
                    self.significant_components[component] = self.normalized_per_component_states[component]

            if to_save:
                for i in self.significant_components:
                    f = self.path_to_results+'/significant_components_' + str(
                        i) + '.npy'
                    np.save(f, self.significant_components[i])




    def run_causal_inference(self, no_combination=None):
        '''
        The function runs the causal analysis on the extracted components. If the n.o. components of choice exceeds 3,
         all of the possible combinations by no_combination variable (recomennded = 3) are checked from the list.
        :param no_combination: int
        :return:
        '''
        if no_combination:
            possible_combinations = [p for p in multiset_permutations([i for i in self.significant_components.keys()],
                                                                      no_combination)]
        else:
            possible_combinations = [p for p in multiset_permutations([i for i in self.significant_components.keys()])]

        results = {}
        taus = [i for i in range(self.window_to_detect)]

        pval_min_threshold, pval_max_threshold = 0.08, 0.08

        frequent_components = {}
        for c in self.significant_components:
            frequent_components[c] = 0

        for combination in possible_combinations:
            print('combination', combination)
            key_combination = tuple(combination)

            results[key_combination] = {}

            for tau in tqdm(taus):

                for t0 in range(self.window_to_detect * 2 - tau * 2):

                    x = self.significant_components[key_combination[0]][:, t0]
                    y = self.significant_components[key_combination[1]][:, t0 + tau * 2]
                    z = self.significant_components[key_combination[2]][:, t0 + tau]

                    data1 = pd.DataFrame({'x': x, 'y': y, 'covar': z})

                    pc = partial_corr(data1, 'x', 'y', 'covar')

                    correlation = corr(data1['x'], data1['y']).round(3)

                    if pc['p-val'][0] > pval_max_threshold:
                        if correlation['p-val'][0] < pval_min_threshold:
                            frequent_components[key_combination[0]] +=1

                            if tau in results[key_combination].keys():
                                results[key_combination][tau].append(t0)
                            else:

                                results[key_combination][tau] = [t0]

        self.per_combination_frequency = results
        self.frequent_components = frequent_components

    def get_center_channels(self):

        component_to_freq = [(c, self.frequent_components[c]) for c in self.frequent_components]
        component_to_freq.sort(key = lambda x:x[1], reverse=True)
        chosen = component_to_freq[:3]
        chosen_indxs = [int(i[0][3:]) for i in chosen]

        per_channel_freq = {}

        for component in chosen_indxs:
            component_mapping = self.topomap[:, component]
            top_5_indices = np.argpartition(component_mapping, -5)[-5:]
            for i in top_5_indices:
                ch = self.ch_names[i]
                if ch in per_channel_freq.keys():
                    per_channel_freq[ch] += 1
                else:
                    per_channel_freq[ch] = 1

        channel_to_freq = [(c, per_channel_freq[c]) for c in per_channel_freq]
        channel_to_freq.sort(key = lambda x:x[1], reverse=True)

        if len(channel_to_freq) > 5:
            chosen_channels = channel_to_freq[:5]
            chosen_channels = [i[0] for i in chosen_channels]
        else:
            chosen_channels = channel_to_freq
            chosen_channels = [i[0] for i in chosen_channels]

        return chosen_channels


    def visualize_per_Tau_frequency(self, tau, savefigures = False):
        '''
        The function plots tick-plot for every causal event detected relative to an IED time for every possible combination of channels.
        :param tau: int
        :return:
        '''
        fig = plt.figure(figsize=(10, 10))
        timescale = list([i for i in range(-self.window_to_detect, self.window_to_detect)])
        plt.title("τ = " + str(tau))
        fixed_tau = tau
        Y_Tick_List = []
        Y_Tick_Label_List = []

        for indx, schema in enumerate(self.per_combination_frequency):

            plt.plot(timescale, list([indx * 10 for k in range(2*self.window_to_detect)]), label=schema[0] + ' -> ' + schema[2])
            Y_Tick_List.append(indx * 10)
            Y_Tick_Label_List.append(schema[0] + ' -> ' + schema[2])

            if fixed_tau in self.per_combination_frequency[schema].keys():
                all_markings = set([k - self.window_to_detect for k in self.per_combination_frequency[schema][fixed_tau]])
                plt.scatter(list(all_markings), list([indx * 10 for k in range(len(all_markings))]), marker="|")
            else:
                continue

        plt.yticks(ticks=Y_Tick_List, labels=Y_Tick_Label_List, rotation=25, fontsize=8)
        plt.axvline(x=0)

        if savefigures:
            plt.savefig(self.path_to_visual_outputs + '/per_tau_results/tau_'+str(tau)+'.jpg')
        plt.close()
