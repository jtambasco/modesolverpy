import numpy as np

def _make_gaussian(x_pts, y_pts, mfd, x_offset=0, y_offset=0):
    x0 = (x_pts[-1]+x_pts[0])/2 + x_offset
    y0 = (y_pts[-1]+y_pts[0])/2 + y_offset
    xx, yy = np.meshgrid(x_pts, y_pts)

    sigma = mfd * 0.707 / 2.355
    sigma_x = sigma
    sigma_y = sigma

    gaus_2d = np.exp(-((xx-x0)**2/(2*sigma_x**2)+
                       (yy-y0)**2/(2*sigma_y**2)))
    gaus_2d /= np.sum(gaus_2d)

    return gaus_2d

def _overlap(mode, gaussian):
    mode_1 = mode
    mode_2 = np.sqrt(gaussian) # square-root for E-field (not power)
    eta = np.abs(np.sum(np.conj(mode_1)*mode_2))**2 / \
        (np.sum(np.abs(mode_1)**2) * np.sum(np.abs(mode_2)**2))
    return eta

def reflection(n1, n2):
    '''
    Calculate the power reflection at the interface
    of two refractive index materials.

    Args:
        n1 (float): Refractive index of material 1.
        n2 (float): Refractive index of material 2.

    Returns:
        float: The percentage of reflected power.
    '''
    r = abs((n1-n2) / (n1+n2))**2
    return r

def transmission(n1, n2):
    '''
    Calculate the power transmission at the interface
    of two refractive index materials.

    Args:
        n1 (float): Refractive index of material 1.
        n2 (float): Refractive index of material 2.

    Returns:
        float: The percentage of transmitted power.
    '''
    return 1-reflection(n1, n2)

def coupling_efficiency(mode_solver, fibre_mfd,
                        fibre_offset_x=0, fibre_offset_y=0,
                        n_eff_fibre=1.441):
    '''
    Finds the coupling efficiency between a solved
    fundamental mode and a fibre of given MFD.

    Args:
        mode_solver (_ModeSolver): Mode solver that
            has found a fundamental mode.
        fibre_mfd (float): The mode-field diameter
            (MFD) of the fibre.
        fibre_offset_x (float): Offset the fibre
            from the centre position of the window
            in x. Default is 0 (no offset).
        fibre_offset_y (float): Offset the fibre
            from the centre position of the window
            in y. Default is 0 (no offset).
        n_eff_fibre (float): The effective index
            of the fibre mode.  Default is 1.441.

    Returns:
        float: The power coupling efficiency.
    '''
    etas = []

    gaus = _make_gaussian(mode_solver._structure.xc, mode_solver._structure.yc,
                          fibre_mfd, fibre_offset_x, fibre_offset_y)

    for mode, n_eff in zip(mode_solver.modes, mode_solver.n_effs):
        o = abs(_overlap(mode, gaus))
        t = abs(transmission(n_eff, n_eff_fibre))
        eta = o * t
        etas.append(eta)

    return etas
