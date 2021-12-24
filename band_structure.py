"""
Program tests array of photon frequencies and produces text file with boolean information about whether a
frequency is allowed for forbidden for a given DCE scattering sideband.

This is the very first computational step in calculating the output DCE radiation of a SQUID-CPW periodic lattice.
Allowed frequencies in main and side bands determine the nature of scattering processes in the lattice.

"""

import numpy as np
import j_constants
from input_frequencies import input_frequencies

j = np.complex(0, 1)  # Define complex number j for convenience


def allowed_frequencies(q_input_freq, N_bins, N_sidebands):
    """
        Performs test to determine allowed and forbidden frequencies for bloch modes in the KP model.

        :param q_input_freq = 2D array with input frequencies
        :param n_sidebands = Number of sidebands used in numerical calculation.
        :param n_bins = Number of bins to break up the frequency range into.
        :return allowed_freq = 2D array with bool values, holds information about whether corresponding element in
                                                input frequency kk array is allowed or forbidden.
        """

    kk = np.abs(q_input_freq)
    condition = np.abs((np.cos(kk)) + ((epsilon / (2 * kk)) * np.sin(kk)))
    allowed_freqs = np.empty((N_bins - 1, 2 * N_sidebands + 1))

    for k in range(0, N_bins - 1):

        # Check center band
        # Note: the index for    the main band is the same as N_sidebands since indices run from zero
        if condition[k, N_sidebands] < 1:
            # If center band IS allowed, make center value True
            allowed_freqs[k, N_sidebands] = True

            # Do forward and backward checks:

            # Check forward:
            for band in range(N_sidebands + 1, 2 * N_sidebands + 1):
                if condition[k, band] < 1:
                    # If band is allowed, make test True
                    allowed_freqs[k, band] = True
                else:
                    # If band is not allowed, make all FORWARD bands False, and don't check following bands (break loop)
                    allowed_freqs[k, band:] = False
                    # print("Forbidden band:", band - N_sidebands, " for bin:", k)
                    break

            # Check backward:
            # Go from center_band-1 to 0 (doesn't include end) in -1 increments
            for band in range(N_sidebands - 1, -1, -1):
                if condition[k, band] < 1:
                    # If band IS allowed, make test True
                    allowed_freqs[k, band] = True
                else:
                    # If band IS NOT allowed, make all BACKWARD bands False, and don't check previous bands (break loop)
                    allowed_freqs[k, band::-1] = False
                    # print("Forbiden band:", band - N_sidebands, " for bin:", k)
                    break
        else:
            # If center band IS NOT allowed, make test False for every band
            allowed_freqs[k, :] = False
        # print("Forbiden band: 0  for bin:", k)

    return allowed_freqs


# Calculation parameters

N_sidebands = 5  # number of sidebands used for calculation
N_bins = 1000  # number of bins dividing (0, 2*N_Omega + 1)

# Physical parameters

ell = 14*j_constants.L0_eff  # lattice separation (in meters)
epsilon = ell/j_constants.L0_eff  # unitless lattice parameter
physical_freq = 18.6  # physical drive frequency of SQUIDs (in GHz)
QQ = np.around(2*np.pi*physical_freq*(10**9)*(ell/j_constants.vcpw), decimals=8)  # unitless drive frequency parameter

"""
----------------------------------------------------------------------------------------------------------------------
MAIN BODY
----------------------------------------------------------------------------------------------------------------------
"""

freq_range = np.array((0, 1))  # Bounds for the input frequencies array, divided by drive frequency.

input_freqs, input_shape = input_frequencies(False, freq_range, N_bins, N_sidebands, QQ)

# Creates grid of input frequencies

allowed_freqs = allowed_frequencies(input_freqs, N_bins, N_sidebands)

np.savetxt("band_structure.csv", allowed_freqs, fmt='% 1d')
