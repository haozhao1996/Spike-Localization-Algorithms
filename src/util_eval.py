import MEArec as mr
import numpy as np
import scipy.optimize
import os
import sys
import re
import ast
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict

import spikeinterface as si
import spikeinterface.core as sc
import spikeinterface.extractors as se
import spikeinterface.sorters as ss
import spikeinterface.preprocessing as spre
import spikeinterface.postprocessing as spost
import spikeinterface.widgets as sw
import spikeinterface.comparison as scomp
import spikeinterface.curation as scu

from util_loc import *

def get_recording_noise(recording, channel_ids, savepath):
    """
    Replaces selected channels in a RecordingExtractor object with Gaussian noise.

    Parameters
    --------------------------
    recording: RecordingExtractor object without noise
    channel_ids: list of int representing dead channels
    
    Returns
    --------------------------
    None (saves the noisy recording in the output folder)
    """
    sampling_frequency = recording.get_sampling_frequency()
    traces = recording.get_traces().copy()
    
    for channel_id in channel_ids:
        traces[:, channel_id] = np.random.normal(loc=0., scale=50., size=traces.shape[0])

    if not os.path.isfile(f'{savepath}/binary.json'):
        recording_numpy = se.NumpyRecording(traces_list=traces, sampling_frequency=sampling_frequency)
        
        # Copy all properties from the original recording
        for prop_name in recording.get_property_keys():
            recording_numpy.set_property(prop_name, recording.get_property(prop_name))
        
        # Set probe (geometry and contact information)
        if recording.get_probe() is not None:
            recording_numpy.set_probe(recording.get_probe())
        
        # Apply preprocessing
        recording_numpy = spre.bandpass_filter(recording_numpy, freq_min=300, freq_max=3000)
        recording_numpy = spre.common_reference(recording_numpy, operator='median', reference='global')
        
        bad_channel_ids = spre.detect_bad_channels(recording=recording_numpy)
        recording_numpy = recording_numpy.remove_channels(remove_channel_ids=bad_channel_ids)

        recording_numpy.save(folder=savepath)
        
    return

def get_unit_loc_est(method, wes):
    """
    Compute the estimated template locations for each day

    Parameters
    --------------------------
    method: str, method to compute unit locations
    wes: list of WaveformExtractor objects

    Returns
    --------------------------
    unit_loc_est: list of numpy arrays, estimated template locations for each day
    """

    unit_loc_est = []
    
    if method == 'dipolar_triangulation':
        for we in wes:
            unit_loc_est.append(compute_dipolar_triangulation(we))
    elif method == 'epm':
        for we in wes:
            unit_loc_est.append(compute_epm_templates(we, optimizer='least_square', radius_um=75, max_distance_um=100, return_alpha = False, enforce_decrease = False, feature='energy'))
    elif method == 'monopolar_triangulation':
        for we in wes:
            unit_loc_est.append(spost.compute_unit_locations(we, method=method, optimizer='least_square'))
    else:        
        for we in wes:
            unit_loc_est.append(spost.compute_unit_locations(we, method=method))

    return unit_loc_est

def get_unit_loc_grid(method, wes, grid_params, days):
    """
    For grid search. Only calculate results for D0, D8, D16, and D24 to save time.
    """

    # Define which days to process
    target_days = ['D0', 'D8', 'D16', 'D24']    
    unit_loc_est = []
    
    for we_i, we in enumerate(wes):
        # Check if current day should be processed
        current_day = days[we_i]
        
        if current_day in target_days:
            # Process this day
            if method == 'monopolar_triangulation':
                unit_locs = spost.compute_unit_locations(we, method=method, optimizer='least_square', radius_um=grid_params['radius_um'], max_distance_um=grid_params['max_distance_um'])
            elif method == 'grid_convolution':
                unit_locs = spost.compute_unit_locations(we, method=method, radius_um=grid_params['radius_um'], percentile=grid_params['percentile'])
            else:  # center_of_mass or other methods
                unit_locs = spost.compute_unit_locations(we, method=method)
            
            unit_loc_est.append(unit_locs)
    
    return unit_loc_est

def get_loc_est_spikes(method, wes):
    """
    Compute the estimated spike locations for each day

    Parameters
    --------------------------
    method: str, method to compute unit locations
    wes: list of WaveformExtractor objects

    Returns
    --------------------------
    unit_loc_est: list of [_], estimated spike locations for each day
    """

    spike_loc_est = []
    
    if method == 'dipolar_triangulation':
        for we in wes:
            spike_loc_est.append(compute_dipolar_triangulation(we))
    elif method == 'epm':
        for we in wes:
            spike_loc_est.append(compute_epm_spikes(we, optimizer='least_square', radius_um=75, max_distance_um=100, return_alpha = False, enforce_decrease = False, feature='energy'))
    elif method == 'monopolar_triangulation':
        for we in wes:
            spike_loc_est.append(spost.compute_spike_locations(we, method=method, outputs='by_unit', max_distance_um=10000))
    else:        
        for we in wes:
            spike_loc_est.append(spost.compute_spike_locations(we, method=method, outputs='by_unit'))

    return spike_loc_est

def get_loss_rmse_gt(unit_loc_true, unit_loc_est):
    """
    Compute the RMSE between true and estimated unit locations using GT sorting results

    Parameters
    --------------------------
    unit_loc_true: numpy array, true unit locations
    unit_loc_est: numpy array, estimated unit locations

    Returns
    --------------------------
    loss: float, RMSE between true and estimated unit locations
    """

    diff = unit_loc_true[:, 1:] - unit_loc_est[:, :2]
    loss = np.sqrt(np.sum(diff**2) / diff.shape[0])

    return loss

def get_loss_mae_gt(unit_loc_true, unit_loc_est):
    """
    Compute the median absolute error between true and estimated unit locations using GT sorting results.
    Cap the error per unit at 500 microns, so one unit does not dominate the loss.

    Parameters
    --------------------------
    unit_loc_true: numpy array, true unit locations
    unit_loc_est: numpy array, estimated unit locations

    Returns
    --------------------------
    loss: float, Median AE between true and estimated unit locations
    """

    diff = unit_loc_true[:, 1:] - unit_loc_est[:, :2]
    distances = np.sqrt(np.sum(diff**2, axis=1))
    for i in range(len(distances)):
        distances[i] = np.minimum(distances[i], 500)
    loss = np.median(distances)

    return loss

def get_loss_median_gt(unit_loc_true_list, unit_loc_est_list):
    """
    Compute the median absolute error across multiple datasets (seeds).

    Parameters
    --------------------------
    unit_loc_true_list: list of numpy arrays, true unit locations for each seed
    unit_loc_est_list: list of numpy arrays, estimated unit locations for each seed

    Returns
    --------------------------
    loss: float, Median AE across all seeds
    """
    
    all_distances = []
    
    for unit_loc_true, unit_loc_est in zip(unit_loc_true_list, unit_loc_est_list):
        diff = unit_loc_true[:, 1:] - unit_loc_est[:, :2]
        distances = np.sqrt(np.sum(diff**2, axis=1))
        all_distances.extend(distances)
    
    # Now compute median across all seeds combined
    loss = np.median(all_distances)
    
    return loss

def get_loss_RMSE_spikes_2(unit_loc_true, unit_loc_est):
    """
    TBU
    """
    
    loss = []
    for spike_loc_est in unit_loc_est:
        spike_loc_est = np.array(spike_loc_est.tolist())[:2]
        diff = unit_loc_true[1:] - spike_loc_est
        loss.append(np.sum(diff**2))

    unit_SE = np.sum(loss)
    unit_spikes = len(loss)

    return unit_SE, unit_spikes

def get_loss_MAE_spikes_2(unit_loc_true, unit_loc_est):
    """
    Mean absolute error
    Cap the error per unit at 500 microns, so one unit does not dominate the loss.
    """
    
    loss = []
    for spike_loc_est in unit_loc_est:
        spike_loc_est = np.array(spike_loc_est.tolist())[:2]
        diff = unit_loc_true[1:] - spike_loc_est
        distance = np.sqrt(np.sum(diff**2))
        distance = np.minimum(distance, 500)  # Cap at 500 microns
        loss.append(distance)

    unit_SE = np.sum(loss)
    unit_spikes = len(loss)

    return unit_SE, unit_spikes

def get_loss_median_spikes(unit_loc_true_list, unit_loc_est_list):
    """
    Compute the median absolute error across multiple datasets (seeds) for spikes.
    """
    
    all_distances = []
    
    for unit_loc_true, unit_loc_est in zip(unit_loc_true_list, unit_loc_est_list):
        for spike_loc_est in unit_loc_est:
            spike_loc_est = np.array(spike_loc_est.tolist())[:2]
            diff = unit_loc_true[1:] - spike_loc_est
            distance = np.sqrt(np.sum(diff**2))
            all_distances.append(distance)
    
    # Now compute median across all spikes from all seeds combined
    loss = np.median(all_distances)
    
    return loss

def get_loss_RMSE_spikes_ind(unit_loc_true, unit_loc_est):
    """
    TBU
    """
    
    loss = []
    for spike_loc_est in unit_loc_est:
        spike_loc_est = np.array(spike_loc_est.tolist())[:2]
        diff = unit_loc_true[1:] - spike_loc_est
        loss.append(np.sqrt(np.sum(diff**2)))

    return loss

def get_loss_acc_spikes_2(unit_loc_true, unit_loc_est, correct_radius):
    """
    Helper function for get_loss_acc_spikes()
    """
    
    n_right = 0
    n_wrong = 0
    for spike_loc_est in unit_loc_est:
        spike_loc_est = np.array(spike_loc_est.tolist())[:2]
        if np.linalg.norm(spike_loc_est - unit_loc_true[1:]) < correct_radius:
            n_right += 1
        else:
            n_wrong += 1

    return n_right, n_wrong

def get_loss_acc_spikes(method, days, correct_radius, unit_loc_true, wes, plot=True):
    """
    Computes the accuracy of a localization method on spikes (not templates)    
    """

    loss = []
    for we_i, we in enumerate(wes):
        print(method, we_i)

        n_right = 0
        n_wrong = 0
        counter = 0
        while n_right < 1 and counter < 10:
            counter += 1
            
            # Compute spike locations
            if method == 'epm':
                unit_loc_est = compute_epm_spikes(we)
            else:
                unit_loc_est = spost.compute_spike_locations(we, method=method, outputs='by_unit')
            
            # Calculate loss
            for unit_idx, unit_id in enumerate(we.unit_ids):
                unit_n_right, unit_n_wrong = get_loss_acc_spikes_2(unit_loc_true[unit_idx], unit_loc_est[0][unit_id], correct_radius)
                n_right += unit_n_right
                n_wrong += unit_n_wrong

        loss.append(n_right / (n_right + n_wrong))

    # Create plot
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))    
        ax.scatter(days, loss)
        ax.title.set_text(f'Spike Accuracy - {method}')
        ax.set_xlabel('Day Number')
        ax.set_ylabel('Accuracy (% Predictions within Tolerance)')

    return loss

def get_loss_acc_templates(method, days, correct_radius, unit_loc_true, wes, plot=True):
    """
    Computes the accuracy of a localization method on templates (not spikes)
    """

    unit_loc_est = get_unit_loc_est(method, wes)

    loss = []
    for we_i, we in enumerate(wes):
        # Calculate loss
        n_right = 0
        n_wrong = 0
        for unit_idx, unit_id in enumerate(we.unit_ids):
            if np.linalg.norm(unit_loc_est[we_i][unit_idx, :2] - unit_loc_true[unit_idx, 1:]) < correct_radius:
                n_right += 1
            else:
                n_wrong += 1
        loss.append(n_right / (n_right + n_wrong))

    # Create plot
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))    
        ax.scatter(days, loss)
        ax.title.set_text(f'Template Accuracy - {method}')
        ax.set_xlabel('Day Number')
        ax.set_ylabel('Accuracy (% Predictions within Tolerance)')

    return loss

def compare_loss_acc(losses, methods, days):
    """
    Creates line plot comparing accuracy (on spikes or templates) across days for multiple methods
    """

    fig, ax = plt.subplots(1, 1, figsize=(12, 10))    
    for method in methods:

        loss = losses[method]

        # Create plot
        ax.plot(days, loss, label=f'{method}')
        ax.title.set_text('Accuracy (% Predictions within Tolerance)')
        ax.set_xlabel('Day Number')
        ax.set_ylabel('Accuracy (% Predictions within Tolerance)')
        ax.legend()

def compare_spike_template(method, unit_loc_true, electrode_loc, dead_indices, wes):
    """
    Creates scatter plot comparing spike (density) and templates (point) location estimates against GT
    """
    range_x = electrode_loc[:, 0].max() - electrode_loc[:, 0].min()
    range_y = electrode_loc[:, 1].max() - electrode_loc[:, 1].min()
    range_A = min(range_x, range_y)

    col_num = 10
    row_num = max(int(np.ceil(len(wes)/col_num)),2)
    fig, axs = plt.subplots(row_num, col_num, figsize=(5*col_num, 5*row_num))  

    for we_i, we in enumerate(wes):
        
        loc_est_spike = spost.compute_spike_locations(we, ms_before=0.1, ms_after=0.1, method=method, outputs='by_unit')
        loc_est_unit = spost.compute_unit_locations(we, method=method)
        
        dead_electrodes = dead_indices[we_i]
        live_electrodes = np.delete(np.arange(electrode_loc.shape[0]), dead_electrodes)

        # Create plot
        row_gt = we_i // col_num
        col_gt = we_i % col_num

        for unit_idx, unit_id in enumerate(we.unit_ids):
            unit_loc_est_temp = loc_est_spike[0][unit_id]
            unit_loc_est_temp = np.array(unit_loc_est_temp.tolist())
            axs[row_gt, col_gt].scatter(unit_loc_est_temp[:, 0], unit_loc_est_temp[:, 1], color=plt.get_cmap('terrain')(unit_idx*20), alpha=0.01, label='Estimated Spikes')

        axs[row_gt, col_gt].scatter(loc_est_unit[:, 0], loc_est_unit[:, 1], color='red', label='Estimated Templates')
        axs[row_gt, col_gt].scatter(unit_loc_true[:, 1], unit_loc_true[:, 2], color='green', label='True Neurons')
        axs[row_gt, col_gt].scatter(electrode_loc[live_electrodes, 0], electrode_loc[live_electrodes, 1], color=[0, 0, 0], alpha=0.1 , label='Electrodes')
        axs[row_gt, col_gt].set_xlim([electrode_loc[:, 0].min() - range_A*0.5, electrode_loc[:, 0].max() + range_A*0.5])
        axs[row_gt, col_gt].set_ylim([electrode_loc[:, 1].min() - range_A*0.5, electrode_loc[:, 1].max() + range_A*0.5])
        axs[row_gt, col_gt].get_xaxis().set_visible(False)
        axs[row_gt, col_gt].get_yaxis().set_visible(False)
        axs[row_gt, col_gt].set_aspect('equal')

def plot_loss_1(method, days, unit_loc_true, sorting_gt, sorting_ex, wes_gt, wes_ex):
    """
    Plots template rmse loss using get_match_all()
    """

    # Compute unit locations
    unit_loc_est_gt = get_unit_loc_est(method, wes_gt)
    unit_loc_est_ex = get_unit_loc_est(method, wes_ex)

    match_ex_gt = get_match_all(sorting_gt, sorting_ex)

    # Compute losses
    loss_gt = []
    loss_ex = []
    for day_i, day in enumerate(days):
        loss_gt.append(get_loss_rmse_gt(unit_loc_true, unit_loc_est_gt[day_i]))
        loss_ex.append(get_loss_rmse_ex(unit_loc_true, unit_loc_est_ex[day_i], match_ex_gt))
        
    # Create plot
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))    
    ax.scatter(days, loss_gt, label='GT Sorting - RMSE')
    ax.scatter(days, loss_ex, label='EX Sorting - RMSE')
    ax.title.set_text(f'Root Mean Squared Error - {method}')
    ax.set_xlabel('Day Number')
    ax.set_ylabel('RMSE')
    ax.legend()

def plot_loss_comp(methods, days, unit_loc_true, sortings_gt, sorting_gt, sortings_ex, sorting_ex, wes_gt, wes_ex):
    """
    Plot loss using get_match_all() for multiple methods
    """

    fig, ax = plt.subplots(1, 1, figsize=(12, 10))    
    for method in methods:

        # Compute unit locations
        unit_loc_est_gt = get_unit_loc_est(method, wes_gt)
        unit_loc_est_ex = get_unit_loc_est(method, wes_ex)

        match_ex_gt = get_match_all(sorting_gt, sorting_ex)

        # Compute losses
        loss_gt = []
        loss_ex = []
        for day_i, day in enumerate(days):
            loss_gt.append(get_loss_rmse_gt(unit_loc_true, unit_loc_est_gt[day_i]))
            loss_ex.append(get_loss_rmse_ex(unit_loc_true, unit_loc_est_ex[day_i], match_ex_gt))
            
        # Create plot
        ax.plot(days, loss_gt, label=f'GT Sorting - {method}')
        ax.plot(days, loss_ex, label=f'EX Sorting - {method}')
        ax.title.set_text(f'Root Mean Squared Error')
        ax.set_xlabel('Day Number')
        ax.set_ylabel('RMSE')
        ax.legend()

def remove_L_suffix(s):
    """
    Utility function for reading meta csv.
    """
    return re.sub(r'(\d+)L', r'\1', s)  

def read_csv_1(filepath):
    """
    Read the meta csv file and return a dictionary.
    """
    # Check the DataFrame structure
    df = pd.read_csv(filepath)
    df.head()

    # Initialize the defaultdict
    result_dict = defaultdict(list)

    # Assuming the DataFrame contains one row with two entries
    if not df.empty:
        key_1_str = remove_L_suffix(df.iloc[0, 0])
        value_1_str = remove_L_suffix(df.iloc[0, 1])
        
        result_dict['npx'] = ast.literal_eval(key_1_str)
        result_dict['patch'] = ast.literal_eval(value_1_str)

    # Display the resulting dictionary
    return result_dict

def detect_peak_on_patch_sig(patch_sig, sample_rate):
    """
    Credit to Samuel Garcia: https://spikeinterface.github.io/blog/marques-smith-neuropixel-384ch-paired-recording/
    """
    # Filter because some traces have drift
    sos = scipy.signal.iirfilter(5, 200./sample_rate*2, analog=False, btype = 'highpass', ftype = 'butter', output = 'sos')
    patch_sig_f = scipy.signal.sosfiltfilt(sos, patch_sig, axis=0)
    
    med = np.median(patch_sig_f)
    mad = np.median(np.abs(patch_sig_f-med))*1.4826 # median absolute deviation
    thresh = med - 12 * mad # 12 MADs below median
    
    # 1 ms refractory period
    d = int(sample_rate * 0.001)
    spike_indexes, prop = scipy.signal.find_peaks(-patch_sig_f, height=-thresh, distance=d)

    return spike_indexes