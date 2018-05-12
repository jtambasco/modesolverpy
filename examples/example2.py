import modesolverpy.mode_solver as ms
import modesolverpy.structure as st
import opticalmaterialspy as mat
import numpy as np

wl = 1.55
x_step = 0.06
y_step = 0.06
wg_height = 0.8
wg_width = 1.8
sub_height = 1.0
sub_width = 4.
clad_height = 1.0
film_thickness = 1.2
angle = 60.

def struct_func(n_sub, n_wg, n_clad):
    return st.RidgeWaveguide(wl, x_step, y_step, wg_height, wg_width,
                             sub_height, sub_width, clad_height,
                             n_sub, n_wg, angle, n_clad, film_thickness)

n_sub = mat.SiO2().n(wl)
n_wg_xx = mat.Ktp('x').n(wl)
n_wg_yy = mat.Ktp('y').n(wl)
n_wg_zz = mat.Ktp('z').n(wl)
n_clad = mat.Air().n()

struct_xx = struct_func(n_sub, n_wg_xx, n_clad)
struct_yy = struct_func(n_sub, n_wg_yy, n_clad)
struct_zz = struct_func(n_sub, n_wg_zz, n_clad)

struct_ani = st.StructureAni(struct_xx, struct_yy, struct_zz)
struct_ani.write_to_file()

solver = ms.ModeSolverFullyVectorial(8)
solver.solve(struct_ani)
solver.write_modes_to_file()

solver.solve_ng(struct_ani, 0.01)

solver.solve_sweep_wavelength(struct_ani, np.linspace(1.501, 1.60, 21))
