## Utilities for Microbiome project

import pickle
import numpy as np
import pandas as pd

def load_pickle(mouse_1='mouse_set_1_data.pkl', mouse_2='mouse_set_2_data.pkl'):
    """ simply load the pickle files """
    with open(mouse_1, 'rb') as f:
        mouse_1 = pickle.load(f)

    with open(mouse_2, 'rb') as f:
        mouse_2 = pickle.load(f)

    return mouse_1, mouse_2

mouse_1, mouse_2 = load_pickle()

def get_mouse_data(mouse_num, column):
    """ helper function to get mouse data """
    if mouse_num <= 5:
        # dict 1
        if column in['reads','qpcr']:
            return mouse_1[column][str(mouse_num)]
        else:
            return mouse_1[column]
    else:
        # dict 2
        if column in['reads','qpcr']:
            return mouse_2[column][str(mouse_num)]
        else:
            return mouse_2[column]

def preprocess_mouse(mouse_reads_df, mouse_qpcr):
    """ converts relative abundances to absolute abundances by first
        converting to proportions of bugs and then multiplying by
        QPCR (mean and std)
    """
    mouse_reads = mouse_reads_df.values
    mouse_qpcr = mouse_qpcr.values

    mouse_read_totals = np.sum(mouse_reads, axis=0)
    mouse_read_percent = mouse_reads * 1/mouse_read_totals

    mouse_read_check = np.sum(mouse_read_percent, axis=0)
    assert np.sum(mouse_read_check) == mouse_reads.shape[1], 'mouse reads proportions didnt sum to one'

    mouse_read_absolute_mean = mouse_read_percent * mouse_qpcr[:,0]
    mouse_read_absolute_std = mouse_read_percent * mouse_qpcr[:,1]

    return mouse_read_absolute_mean, mouse_read_absolute_std, mouse_read_percent

def load_data():
    """
    returns a dictionary of mice and their data

    e.g. mice['2']['reads_abs_mean'] retuns absolute reads mean for mouse 2
    """
    mice = {}
    for mouse in range(2,11):
        mouse_reads = get_mouse_data(mouse, 'reads')
        mouse_qpcr = get_mouse_data(mouse, 'qpcr')
        mouse_times = get_mouse_data(mouse, 'times')
        mouse_otus = get_mouse_data(mouse, 'otu_taxonomy')
        mouse_reads_abs_mu, mouse_reads_abs_std, mouse_reads_percent = preprocess_mouse(mouse_reads, mouse_qpcr)
        # standard scaling mean absolute reads across time
        mouse_reads_standardized = mouse_reads_abs_mu - np.mean(mouse_reads_abs_mu, axis=1).reshape(-1,1)
        mouse_reads_standardized = mouse_reads_standardized / np.std(mouse_reads_abs_mu, axis=1).reshape(-1,1)
        mice[mouse] = dict(otus=np.array(mouse_otus.index),
                          reads=mouse_reads,
                          qpcr=mouse_qpcr,
                          times=np.array(mouse_times),
                          reads_abs_mean=mouse_reads_abs_mu,
                          reads_abs_std=mouse_reads_abs_std,
                          reads_percent=mouse_reads_percent,
                          reads_standardized=mouse_reads_standardized)
    return mice
