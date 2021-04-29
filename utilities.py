import numpy as np
import pywt
import matplotlib.pyplot as plt

def get_coef_count_in_each_dwt_level(p, w, level, mode='symmetric', printon=True):
    """
    Return 1) the number of wavelet coefficients (or features) in each level of the DWT, and 
    2) index+1 of the last coefficient in each level, in a flattened array of coefficients
    """
    coef = pywt.wavedec(np.arange(p), w, level=level, mode=mode) 
    coef_count_per_level = [len(c) for c in coef]
    lastcoef_idx_per_level = list(np.cumsum(np.array(coef_count_per_level)))
    if printon:
        print('# of coefficients in each DWT level:', coef_count_per_level)
        print('Idx+1 of last feature in each level in flattened array of coefficients:', lastcoef_idx_per_level)    
    return coef_count_per_level, lastcoef_idx_per_level


def perform_dwt_on_signal_array(signal_array, wavelet_type='coif4', mode='symmetric', level=None):
    """
    Perform DWT on 2D array of signals. 
    Returns 2D array of flattened DWT coefficients.
    """
    n, p = signal_array.shape
    if level is None:
        level = pywt.dwt_max_level(data_len=p, filter_len=pywt.Wavelet(wavelet_type).dec_len)
    _, coef_idx_bylevel = get_coef_count_in_each_dwt_level(p, wavelet_type, level, mode)
    DWTcoef_bank_list = []
    DWTcoef_array = np.zeros((n, coef_idx_bylevel[-1]))
    
    # iterate through samples 
    for i in range(n):
        x = signal_array[i,:]
        coef = pywt.wavedec(x, wavelet_type, level=level, mode=mode) 
        DWTcoef_bank_list.append(coef)
        # flatten coeffs into vec
        coef_flattened, _ = pywt.coeffs_to_array(coef)
        # update array of DWT features
        DWTcoef_array[i,:] = np.array(coef_flattened)
    return DWTcoef_array


def coef_flattened_to_bank(coef_flattened, p=992, wavelet_type='coif4', mode='symmetric', level=None): 
    """
    Convert a flattened array of DWT coefficients back into a list of arrays, corresponding to the coefficients for each level
    """
    if level is None:
        level = pywt.dwt_max_level(data_len=p, filter_len=pywt.Wavelet(wavelet_type).dec_len)
    _, coef_idx_bylevel = get_coef_count_in_each_dwt_level(p, wavelet_type, level, mode, printon=False)
    coef_idx_bylevel = [0] + coef_idx_bylevel
    coef_bank = []
    for i in range(level+1):
        coef_bank.append(np.array(coef_flattened[coef_idx_bylevel[i]:coef_idx_bylevel[i+1]]))
    return coef_bank


def perform_idwt_on_coef_array(DWTcoef_array, p=992, wavelet_type='coif4', mode='symmetric', level=None):
    """
    Perform inverse DWT to convert an array of flattened DWT coefficients to an array of signals
    """
    reconstructed_signal_array = []
    for i in range(len(DWTcoef_array)):
        coef_flattened = DWTcoef_array[i,:]
        coef_bank = coef_flattened_to_bank(coef_flattened, p, wavelet_type, mode, level=None)
        sig = pywt.waverec(coef_bank, wavelet_type)
        reconstructed_signal_array.append(sig)
    reconstructed_signal_array = np.array(reconstructed_signal_array)
    return reconstructed_signal_array
    
    
def visualize_reconstructed_signal_bylevel(coef_bank, signal=None, wv=None, wavelet_type='coif4', mode='symmetric', plot_approximation_bkgd=True):
    """
    Given DWT coefficient bank (list of array of DWT coefficients for each level), plot 
    """
    lvl = len(coef_bank)-1
    # initialize filtered coefficient bank
    coef_bkgd = []
    coef_filt = []
    for coefs_lvl in coef_bank:
        coef_bkgd.append(np.zeros(len(coefs_lvl),))
    if plot_approximation_bkgd:
        coef_filt = coef_bkgd.copy()
        coef_filt[0] = coef_bank[0]
    
    fig, ax = plt.subplots(int(np.ceil((lvl+1)/3)), 3, sharex='col', sharey='row', figsize=(16,10))
    color = ['k','m','b','g','orange','r']
    
    # plot the various levels of coefficients 
    for l in range(lvl+1):
        # get plot indices
        i = int(np.floor(l/3))
        j = l%3
        if not plot_approximation_bkgd:
            coef_filt = coef_bkgd.copy()
        coef_filt[l] = coef_bank[l]
        if l==0:
            plot_title = 'Signal reconstruction using cA' + str(lvl) + ' coefficients'
        else:
            if plot_approximation_bkgd:
                plot_title = 'Signal reconstruction using cD' + str(lvl-l+1) + ' + cA5 coefficients'
            else: 
                plot_title = 'Signal reconstruction using cD' + str(lvl-l+1) + ' coefficients'
        sig_rec_partial = pywt.waverec(coef_filt, wavelet_type, mode=mode, axis=-1)
        if wv is None: 
            wv = np.arange(len(sig_rec_partial))
        if signal is not None:
            ax[i,j].plot(wv,signal, linewidth = 0.25, color = 'grey')
        ax[i,j].plot(wv,sig_rec_partial, linewidth=1, color = color[l-1])
        ax[i,j].title.set_text(plot_title)
    plt.show()
    
    
def get_filtered_coef_bank(coef_bank, coef_idx_to_plot):
    """
    Set all DWT coefficients in bank to zero, except for those indices in the flattened coefficient array to be plotted as reconstructed wavelets
    """
    coef_flattened, _ = pywt.coeffs_to_array(coef_bank)
    coef_flattened_filt = np.zeros(len(coef_flattened),)
    # set all coefficients to zero, except those in list of indices to plot
    coef_flattened_filt[np.array(coef_idx_to_plot)] = coef_flattened[np.array(coef_idx_to_plot)] 
    
    # reshape into filtered coefficient bank
    coef_bank_filt = []
    start_idx=0
    for i, coefs_lvl in enumerate(coef_bank):
        num_coefs_lvl = len(coefs_lvl)
        end_idx = start_idx + num_coefs_lvl
        if i < len(coef_bank)-1:
            coef_bank_filt.append(coef_flattened_filt[start_idx:end_idx])
        else:
            coef_bank_filt.append(coef_flattened_filt[start_idx:])
        start_idx = end_idx
    return coef_bank_filt
    
    
def plot_individual_reconstructed_wavelets(coef_bank, signal, wv, coef_idx_to_plot, wavelet_type='coif4', mode='symmetric'): 
    """
    Plot contribution of individual reconstructed wavelets from their respective DWT coefficients, overlaying the original signal
    """
    # get coefficient bank filtered by indices to plot
    coef_bank_filt = get_filtered_coef_bank(coef_bank, coef_idx_to_plot)
    
    # reconstruct wavelet signal based on filtered coefficient bank
    sig_rec_partial = pywt.waverec(coef_bank_filt, wavelet_type, mode=mode, axis=-1)
    
    # plot reconstructed and original signal
    plt.plot(wv, signal)
    plt.plot(wv, sig_rec_partial)
    plt.legend(['original signal', 'filtered/reconstructed signal'])
    plt.show()
    