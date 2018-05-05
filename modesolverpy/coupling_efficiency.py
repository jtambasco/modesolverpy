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
    gaus_sqrt = np.sqrt(gaussian)
    mode_norm = mode / np.sum(mode)
    np.savetxt('gaus.dat', gaussian)
    np.savetxt('gaus_sqrt.dat', gaus_sqrt)
    eta = np.sum(mode_norm*gaus_sqrt)**2 / (np.sum(mode_norm**2) * np.sum(gaus_sqrt**2))
    return eta

def reflection(n1, n2):
    r = abs((n1-n2) / (n1+n2))**2
    return r

def transmission(n1, n2):
    return 1-reflection(n1, n2)

def coupling_efficiency(mode_solver, fibre_mfd,
                        fibre_offset_x=0, fibre_offset_y=0):
    etas = []

    #fibre_mfd /= 1.699

    gaus = _make_gaussian(mode_solver._structure.xc, mode_solver._structure.yc,
                          fibre_mfd, fibre_offset_x, fibre_offset_y)

    for mode, n_eff in zip(mode_solver.modes, mode_solver.n_effs):
        o = abs(_overlap(mode, gaus))
        t = abs(transmission(n_eff, 1.441))
        eta = o * t
        etas.append(eta)

    return etas
