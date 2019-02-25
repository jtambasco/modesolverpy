import modesolverpy.mode_solver as ms
import modesolverpy.structure as st
import opticalmaterialspy as mat
import numpy as np

wl = 1.55
x_step = 0.02
y_step = 0.02
etch_depth = 0.22
wg_widths = np.arange(0.3, 2., 0.02)
sub_height = 1.
sub_width = 4.
clad_height = 1.
film_thickness = 0.22

n_sub = mat.SiO2().n(wl)
n_clad = mat.Air().n(wl)
n_wg = mat.RefractiveIndexWeb(
    'https://refractiveindex.info/?shelf=main&book=Si&page=Li-293K').n(wl)

r = []
for w in wg_widths:
    r.append(
        st.RidgeWaveguide(wl, x_step, y_step, etch_depth, w, sub_height,
                          sub_width, clad_height, n_sub, n_wg, None, n_clad,
                          film_thickness))

r[0].write_to_file('start_n_profile.dat')
r[-1].write_to_file('end_n_profile.dat')

solver = ms.ModeSolverFullyVectorial(6)
solver.solve_sweep_structure(r, wg_widths, x_label='Taper width', fraction_mode_list=[1,2])
solver.write_modes_to_file()
