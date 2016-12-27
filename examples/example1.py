import modesolverpy.mode_solver as ms
import modesolverpy.structure as st

# All units are relative.  [um] were chosen in this case.
x_step = 0.02
y_step = 0.02
wg_height = 0.4
wg_width = 0.5
sub_height = 0.5
sub_width = 2.
clad_height = 0.5
n_sub = 1.4
n_wg = 3.
n_clad = 1.
film_thickness = 0.5
wavelength = 1.55

structure = st.RidgeWaveguide(x_step,
                              y_step,
                              wg_height,
                              wg_width,
                              sub_height,
                              sub_width,
                              clad_height,
                              n_sub,
                              n_wg,
                              n_clad,
                              film_thickness)
structure.write_to_file('example_structure_1.dat')

mode_solver = ms.ModeSolverSemiVectorial(2)
mode_solver.solve(structure, wavelength)
mode_solver.write_modes_to_file('example_modes_1.dat')
