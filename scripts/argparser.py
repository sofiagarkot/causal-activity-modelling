from argparse import ArgumentParser, Namespace
import os

def add_args(parser):
    """
    Adds predict arguments to an ArgumentParser.

    :param parser: An ArgumentParser.
    """
    parser.add_argument('--edf_file', type=str,
                        default='/data/garkots/third_package/pat04_ON_wake/EEG_1277-export.edf',
                        help='Path to edf file with corresponding data for a patient')

    parser.add_argument('--annotations', type=str,
                        default='/data/garkots/third_package/pat04_ON_wake/',
                        help='Path to markings.')

    parser.add_argument('--path_to_visual_outputs',
                        default='/home/garkots/invasive_eeg/analysis/visuals/experiment_window_500/pat4/wake',
                        type=str)

    parser.add_argument('--path_to_results',
                        default='/home/garkots/invasive_eeg/analysis/results/experiment_window_500/pat4/wake',
                        type=str)

    parser.add_argument('--path_to_ica_matrix',
                        default='/home/garkots/invasive_eeg/analysis/matrices_ica/pat4/wake',
                        type=str)

    parser.add_argument('--chosen_components',nargs='+', default=[])

    parser.add_argument('--ignore_channels',
                        default=True,
                        type=bool)

    parser.add_argument('--precomputed_ica',
                        default=False,
                        type=bool)

    parser.add_argument('--control_experiment',
                        default=False,
                        type=bool)

    parser.add_argument('--state', type=str,
                        default='wake',
                        choices = ['wake', 'sleep'])

    parser.add_argument('--window_length', type=int,
                        default=100)

    parser.add_argument('--automated_choice', type=bool,
                        default=False)

    parser.add_argument('--save_visuals', type=bool,
                        default=True)



def check_args(parser):

    if parser.path_to_visual_outputs:
        assert os.path.isdir(parser.path_to_visual_outputs)

        if not os.path.isdir(parser.path_to_visual_outputs+'/all_components'):
            os.mkdir(parser.path_to_visual_outputs+'/all_components')

        if not os.path.isdir(parser.path_to_visual_outputs+'/significant_components'):
            os.mkdir(parser.path_to_visual_outputs+'/significant_components')

        if not os.path.isdir(parser.path_to_visual_outputs+'/per_tau_results'):
            os.mkdir(parser.path_to_visual_outputs+'/per_tau_results')

def parse_args():
    parser = ArgumentParser()
    add_args(parser)
    args = parser.parse_args()
    check_args(args)

    return args


