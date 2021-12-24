"""
Function generates array of frequencies for different scattering sidebands used in numerical DCE calculation.


:param endpoint = bool, whether or not to include endpoint in array (drive frequency)
:param freq_range = 2 element list, bounds of the frequency range examined (in units of drive frequency),
        useful to run diagnostics with limited frequency ranges.
:param n_bins = int, number of bins to break up (0, omega_d) into.
:param n_sidebands = int, number of sidebands used in numerical calculation.
:param omega_d = float, drive frequency (To make q_input unitles, pass QQ as input.)
:return frequency_grid = 2D array, columns for each sideband (ordered n_sidebands to minus n_sidebands) and rows
        for each frequency.
"""

def input_frequencies(endpoint, freq_range, N_bins, N_sidebands, omega_d):
    import numpy as np
    freq_range_low, freq_range_high = freq_range[0], freq_range[1]
    low_index = int(N_bins * freq_range_low) - 1
    high_index = int(N_bins * freq_range_high) - 1

    freq_grid_size = high_index - low_index + 1
    freq_grid_shape = (freq_grid_size, 2*N_sidebands + 1)

    frequency_grid = np.empty(freq_grid_shape, np.double)
    k_space = np.array(range(low_index, high_index + 1), np.double)

    for k in range(0, freq_grid_size):
        for m in range(1, 2*N_sidebands + 2):
            frequency_grid[k, m-1] = (k_space[k]+1) * (omega_d/N_bins) + (N_sidebands + 1 - m) * omega_d

    if not endpoint:
        frequency_grid = frequency_grid[1:-1, :]
        freq_grid_shape = np.shape(frequency_grid)

    return frequency_grid, freq_grid_shape
