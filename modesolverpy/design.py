import numpy as np

def directional_coupler_lc(wavelength_nm, n_eff_1, n_eff_2):
    '''
    Calculates the coherence length (100% power transfer) of a
    directional coupler.
    '''
    wavelength_m = wavelength_nm * 1.e-9
    dn_eff = (n_eff_1 - n_eff_2).real
    lc_m = wavelength_m / (2.*dn_eff)
    lc_um = lc_m * 1.e6
    return lc_um

def grating_coupler_period(wavelength, n_eff, n_clad,
                           incidence_angle_deg, diffration_order=1):

    k0 = 2.*np.pi / wavelength
    beta = n_eff.real * k0
    n_inc = n_clad

    grating_period = (2.*np.pi*diffration_order) \
        / (beta - k0*n_inc*np.sin(np.radians(incidence_angle_deg)))

    return grating_period
