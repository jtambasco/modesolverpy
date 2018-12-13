import modesolverpy.mode_solver as ms
import modesolverpy.structure as st
import modesolverpy.design as de
import opticalmaterialspy as mat
import numpy as np
import tqdm

wls = [1.5, 1.55, 1.6]
x_step = 0.02
y_step = 0.02
etch_depth = 0.22
wg_width = 0.44
sub_height = 0.5
sub_width = 2.
clad_height = 0.5
film_thickness = 0.22
gaps = np.linspace(0.1, 0.5, 11)

for wl in wls:
    lengths = []

    n_sub = mat.SiO2().n(wl)
    n_clad = mat.Air().n(wl)
    n_wg = 3.476

    for gap in tqdm.tqdm(gaps):
        r = st.WgArray(wl, x_step, y_step, etch_depth, [wg_width, wg_width], gap,
                       sub_height, sub_width, clad_height, n_sub, n_wg, None)
        #r.write_to_file()

        solver = ms.ModeSolverFullyVectorial(2)
        solver.solve(r)
        n1 = solver.n_effs_te[0]
        n2 = solver.n_effs_te[1]
        lengths.append(de.directional_coupler_lc(wl*1000, n1, n2)/2)

    filename = 'dc-sweep-%inm-%s-%inm-etch-%i-film.dat' % (wl*1000, 'TE', etch_depth*1000, film_thickness*1000)
    np.savetxt(filename, np.c_[gaps, lengths], delimiter=',', header='Coupling lengths (50\%)')
