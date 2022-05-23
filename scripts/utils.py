import numpy as np
import os
import json
from scipy.stats import zscore
from antropy import perm_entropy
from copy import deepcopy
import matplotlib.pyplot as plt
from tqdm import tqdm

def get_average(arrays):
    # arrays --> list of lists, where every list is a cropped IED on a Signal

    result = []
    stdevs = []

    # iterating over every count
    for enrty_indx in range(len(arrays[0])):

        # we will record average summing the corresponding counts
        to_get_average = []

        # iterating over every component
        for i in range(len(arrays)):
            # get the state
            state = arrays[i]

            to_get_average.append(state[enrty_indx])

        result.append(sum(to_get_average) / len(to_get_average))
        stdevs.append(np.std(to_get_average))

    return np.array(result), np.array(stdevs)


def translate(value, leftMin, leftMax, rightMin, rightMax):

    # Figure out how 'wide' each range is
    leftSpan = leftMax - leftMin
    rightSpan = rightMax - rightMin

    # Convert the left range into a 0-1 range (float)

    valueScaled = float(value - leftMin) / float(leftSpan)

    # Convert the 0-1 range into a value in the right range.
    return rightMin + (valueScaled * rightSpan)



def read_annotations(path):
    f = open(path)

    onsets = {}

    for el in f.readlines()[1:]:
        try:
            time, channel = el.strip().split(',')
        except:
            print(el)
            continue
        channel = str(channel.strip("'"))
        time = float(time)

        if channel not in onsets.keys():
            onsets[channel] = [time]
        else:
            onsets[channel].append(time)

    return onsets

def no_channels_read_annotations(path, fs, window, max_length):

    annotations = {}

    for file in os.listdir(path):
        if file[-4:] == '.evt':

            ch = file[:-4].upper()
            annotations[ch] = []

            f = open(path + file)

            for line in f.readlines()[1:]:
                event = float(line.split()[0]) * 0.000001 * fs
                if int(event) < max_length - window and int(event) - window > 0:
                    annotations[ch].append(int(event))

    # i = 0
    # for c in annotations:
    #     k = annotations[c]
    #
    #     k = sorted(k)
    #     k2 = deepcopy(k)
    #
    #     for indx, el in enumerate(k[:-1]):
    #         if el + int(0.7*window) >= k[indx + 1]:
    #             k2.pop(indx - i)
    #             i += 1
    #
    #     annotations[c] = k2

    return annotations

def deviative_choice(patient):

    per_channel_deviation_names = []
    per_channel_deviation_values  = []
    for channel in tqdm(patient.per_component_state_during_seizures):
        per_channel_deviation_names.append(channel)
        to_shape = np.asarray(patient.per_component_state_during_seizures[channel]).shape[0]*np.asarray(patient.per_component_state_during_seizures[channel]).shape[1]
        concatenated=np.reshape(np.asarray(patient.per_component_state_during_seizures[channel]), to_shape)

        entropy = perm_entropy(concatenated, 3, 1)
        span = np.max(concatenated) - np.min(concatenated)
        combined_val = entropy - span
        per_channel_deviation_values.append(combined_val )

    per_channel_deviation_values = np.asarray(per_channel_deviation_values)
    per_channel_deviation_names = np.asarray(per_channel_deviation_names)

    indx_val = {}
    for i, val in enumerate(per_channel_deviation_values):
        indx_val[val] = i

    print(per_channel_deviation_values)
    chosen = list(per_channel_deviation_names[zscore(per_channel_deviation_values) <= -2])

    while len(chosen) < 3:

        sorted_per_channel_deviation_values = np.sort(per_channel_deviation_values)

        for stdev_val in sorted_per_channel_deviation_values:
            if per_channel_deviation_names[indx_val[stdev_val]] not in chosen:
                chosen.append(per_channel_deviation_names[indx_val[stdev_val]])
                break

    return chosen


def per_tau_events_density_plot(patient):
    '''
    The function visualizes the n.o. causal events for a particular patient related to tau.
    Usage: to determine the tau with the biggest number of causal events.
    :param patient: Patient
    '''

    tau_density = {}
    for i in range(patient.window_to_detect):
        tau_density[i] = 0

    for schema in patient.per_combination_frequency:
        for tau in patient.per_combination_frequency[schema]:
            for event in patient.per_combination_frequency[schema][tau]:
                tau_density[tau] += 1

    plt.figure(figsize=(30, 10))
    plt.plot([i for i in tau_density.keys()], [tau_density[i] for i in tau_density.keys()])
    plt.yticks(fontsize=24)
    plt.xticks(fontsize=24)
    if patient.path_to_visual_outputs:
        plt.savefig(patient.path_to_visual_outputs + '/per_tau_results/tau_'+str(tau)+'.jpg')

    return list([tau_density[i] for i in tau_density.keys()]).index(
        max(list([tau_density[i] for i in tau_density.keys()])))


def get_density_plot(patient, selected_tau=False):

    density = {}
    for i in range(-patient.window_to_detect, patient.window_to_detect):
        density[i] = 0

    total_events = []

    if type(selected_tau) == bool and selected_tau == False:
        for schema in patient.per_combination_frequency:
            for tau in patient.per_combination_frequency[schema]:
                for event in patient.per_combination_frequency[schema][tau]:
                    density[event - patient.window_to_detect] += 1
                    total_events.append(event - patient.window_to_detect)
    else:
        for schema in patient.per_combination_frequency:
            if selected_tau in patient.per_combination_frequency[schema].keys():
                for event in patient.per_combination_frequency[schema][selected_tau]:
                    density[event - patient.window_to_detect] += 1
                    total_events.append(event - patient.window_to_detect)

    plt.figure(figsize=(30, 10))
    plt.hist(total_events, bins=50)
    plt.xlabel('Time relative to seizure', fontsize=24)
    plt.ylabel('№ of causal events detected', fontsize=24)
    plt.yticks(fontsize=24)
    plt.xticks(fontsize=24)

    if type(selected_tau) == bool and selected_tau == False:
        plt.savefig(patient.path_to_results + '/all_taus.jpg')
    else:
        plt.savefig(patient.path_to_results + '/max_tau.jpg')

    return density

def save_to_json(path, dct):

    with open(path, 'w') as f:
        json.dump(dct, f)


def parse_ica_channels():
    channels = []
    to_insert = '000'
    while to_insert!= '---':
        to_insert = str(input('Type component name in format ICA###:'))

        print('Type --- to interrupt...')
        if to_insert != '---':
            channels.append(to_insert)


    return channels