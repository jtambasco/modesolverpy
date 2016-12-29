import numpy as np
import scipy.optimize as sciopt

def gaussian(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def fit_gaussian(x, y, z_2d, save_fits=False):
    z = z_2d

    max_idx = np.unravel_index(z.argmax(), z.shape)
    max_row = max_idx[0] - 1
    max_col = max_idx[1] - 1

    z_max_row = z[max_row, :]
    z_max_col = z[:, max_col]
    A = z[max_row, max_col]

    p_guess_x = (A, x[max_col], 0.1*(x[-1] - x[0]))
    p_guess_y = (A, y[max_row], 0.1*(y[-1] - y[0]))

    coeffs_x, var_matrix_x = sciopt.curve_fit(gaussian, x, z_max_row, p_guess_x)
    coeffs_y, var_matrix_y = sciopt.curve_fit(gaussian, y, z_max_col, p_guess_y)

    c_x = (x[-1]-x[0])*(max_col+1)/x.size + x[0]
    c_y = (y[-1]-y[0])*(y.size-(max_row+1))/y.size + y[0]
    centre = (c_x, c_y)

    sigma = np.array([coeffs_x[2], coeffs_y[2]])
    fwhm = 2.355 * sigma
    sigma_2 = 1.699 * fwhm

    if save_fits:
        with open('x_fit.dat', 'w') as fs:
            for c in np.c_[x, z_max_row, gaussian(x, *coeffs_x)]:
                s = ','.join([str(v) for v in c])
                fs.write(s+'\n')
        with open('y_fit.dat', 'w') as fs:
            for c in np.c_[y, z_max_col, gaussian(y, *coeffs_y)]:
                s = ','.join([str(v) for v in c])
                fs.write(s+'\n')

    return A, centre, sigma_2
