from argparser import parse_args
from analysis import make_patient_analysis

if __name__ == '__main__':
    args = parse_args()
    make_patient_analysis(args)