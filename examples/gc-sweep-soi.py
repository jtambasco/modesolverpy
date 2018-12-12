import modesolverpy.mode_solver as ms
import modesolverpy.structure as st
import modesolverpy.design as de
import opticalmaterialspy as mat
import numpy as np

wls = [1.5, 1.55, 1.6]
x_step = 0.05
y_step = 0.05
etch_depth = 0.07
wg_width = 10
sub_height = 0.5
sub_width = 14.
clad_height = 0.5
film_thickness = 0.22
polarisation = 'TE'
dcs = np.linspace(20, 80, 61) / 100

ed1 = etch_depth
ft1 = film_thickness
ed2 = ft1 - ed1
ft2 = ed2

periods = []
periods.append(dcs)

for wl in wls:
    ngc = []
    for ed, ft in [(ed1, ft1), (ed2, ft2)]:
        def struct_func(n_sub, n_wg, n_clad):
            return st.RidgeWaveguide(wl, x_step, y_step, ed, wg_width,
                                     sub_height, sub_width, clad_height,
                                     n_sub, n_wg, None, n_clad, ft)

        n_sub = mat.SiO2().n(wl)
        n_wg_xx = 3.46
        n_wg_yy = 3.46
        n_wg_zz = 3.46
        n_clad = mat.Air().n()

        struct_xx = struct_func(n_sub, n_wg_xx, n_clad)
        struct_yy = struct_func(n_sub, n_wg_yy, n_clad)
        struct_zz = struct_func(n_sub, n_wg_zz, n_clad)

        struct_ani = st.StructureAni(struct_xx, struct_yy, struct_zz)
        #struct_ani.write_to_file()

        solver = ms.ModeSolverFullyVectorial(4)
        solver.solve(struct_ani)
        #solver.write_modes_to_file()

        if polarisation == 'TE':
            ngc.append(np.round(np.real(solver.n_effs_te), 4)[0])
        elif polarisation == 'TM':
            ngc.append(np.round(np.real(solver.n_effs_tm), 4)[0])

    period = de.grating_coupler_period(wl, dcs*ngc[0]+(1-dcs)*ngc[1], n_clad, 8, 1)
    periods.append(period)

filename = 'dc-sweep-%s-%inm-etch-%i-film.dat' % (polarisation, etch_depth*1000, film_thickness*1000)
np.savetxt(filename, np.array(periods).T, delimiter=',', header=','.join([str(val) for val in wls]))
print(np.c_[periods])
