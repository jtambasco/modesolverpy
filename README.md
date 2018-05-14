# modesolverpy
Photonic mode solver with a nice interface and output.

The documentation for this project can be found [here](http://modesolverpy.rtfd.io).

## Installation
It is recommend to install `modesolverpy` either via:

### Ubuntu/Mint/Debian:
```bash
pip3 install modesolverpy
apt install gnuplot
```

### Arch Linux:
```bash
yaourt -S python-modesolverpy
```

### Dependencies
If installing using the [Arch Linux AUR package](https://aur.archlinux.org/packages/python-modesolverpy/) or `pip`, dependencies will be automatically downloaded and installed, if not, one should ensure the following dependencies are installed:

#### Python

* [setuptools](https://pypi.python.org/pypi/setuptools),
* [numpy](http://www.numpy.org/),
* [scipy](https://www.scipy.org/),
* [tqdm](https://pypi.python.org/pypi/tqdm),
* [gnuplotpy](https://github.com/jtambasco/gnuplotpy), and
* [opticalmaterialspy](https://github.com/jtambasco/opticalmaterialspy).

#### Other

* [Gnuplot](http://www.gnuplot.info/).

## Features
The main reasons to consider this library include:

* semi-vectorial and fully vectorial options,
* simple structure drawing,
* automated data saving and plotting via Gnuplot,
* some limited (at this stage) data processing (finding MFD of fundamental mode), and
* easily extensible library.

## Example 1

### Semi-vectorial mode solving of a ridge waveguide
The following example finds the first two modes of a waveguide with the following, arbitrary, parameters:

* thin-film thickness: 500nm
* waveguide height: 400nm,
* waveguide width: 500nm,
* refractive index of waveguide: 3,
* refractive index of substrate: 1.4,
* refractive index of cladding: 1, and
* wavelength: 1550nm.

#### Python script
```python
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

mode_solver = ms.ModeSolverSemiVectorial(2, semi_vectorial_method='Ey')
mode_solver.solve(structure, wavelength)
mode_solver.write_modes_to_file('example_modes_1.dat')
```

#### Structure
<img src="./examples/example_structure_1.png " width="400">

### Modes
<img src="./examples/example_modes_1_Ey_0.png " width="400"> <img src="./examples/example_modes_1_Ey_1.png " width="400">

## Example 2

### Fully vectorial mode solving of anisotropic material
The following looks at a contrived ridge waveguide in Z-cut KTP.

The simulation outputs:
* 5 plots for each refractive index axis (n_xx, n_xy, n_yx, n_yy and n_zz),
* 48 plots for Ex, Ey, Ez, Hx, Hy and Hz,
* 8 effective index values, one for each mode,
* a wavelength sweep of the waveguide (plotting n_eff vs wavelength for each mode),
* whether a mode is qTE or qTM and the percentage overlap with TE and TM, and
* the group velocity of the mode.

The waveguide parameters are:
* thin-film thickness: 1.2um,
* waveguide height: 800nm,
* waveguide width: 1.2um,
* refractive index of waveguide: used Sellmeier equations to get n_xx, n_yy, n_zz at 1550nm,
* refractive index of substrate: used Sellmeier equation to get SiO2 at 1550nm,
* refractive index of cladding: 1, and
* wavelength: 1550nm.

### Python script
```python
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

solver.solve_ng(struct_ani, 1.55, 0.01)

solver.solve_sweep_wavelength(struct_ani, np.linspace(1.501, 1.60, 21))
```

### Group Velocity
The group velocity at 1550nm for each mode is:
```
# modes_full_vec/ng.dat
# Mode idx, Group index
0,1.776
1,1.799
2,1.826
3,1.847
4,1.841
5,1.882
6,1.872
7,1.871
```

### Structure
<img src="./examples/material_index/material_index_xx.png " width="250"> <img src="./examples/material_index/material_index_xy.png " width="250"> <img src="./examples/material_index/material_index_yx.png " width="250">
<img src="./examples/material_index/material_index_yy.png " width="250"> <img src="./examples/material_index/material_index_zz.png " width="250">

### Modes
Only the first 4 (out of 8) modes are shown, and only the E-fields are shown (not H-fields).  For the rest of the images, look in the example folder or run the script.

A_{x,y,z} give the percentage power of that particular E-field component with respect to the total of all components.

Mode types:
```
# modes_full_vec/mode_info
# Mode idx, Mode type, % in major direction, n_eff
0,qTE,97.39,1.643
1,qTM,92.54,1.640
2,qTE,90.60,1.576
3,qTM,91.41,1.571
4,qTE,89.48,1.497
5,qTM,86.70,1.475
6,qTE,89.47,1.447
7,qTM,68.35,1.437
```

<img src="./examples/modes_full_vec/mode_0/mode_Ex_0.png " width="265"> <img src="./examples/modes_full_vec/mode_0/mode_Ey_0.png " width="265"> <img src="./examples/modes_full_vec/mode_0/mode_Ez_0.png " width="265">
<img src="./examples/modes_full_vec/mode_1/mode_Ex_1.png " width="265"> <img src="./examples/modes_full_vec/mode_1/mode_Ey_1.png " width="265"> <img src="./examples/modes_full_vec/mode_1/mode_Ez_1.png " width="265">
<img src="./examples/modes_full_vec/mode_2/mode_Ex_2.png " width="265"> <img src="./examples/modes_full_vec/mode_2/mode_Ey_2.png " width="265"> <img src="./examples/modes_full_vec/mode_2/mode_Ez_2.png " width="265">
<img src="./examples/modes_full_vec/mode_3/mode_Ex_3.png " width="265"> <img src="./examples/modes_full_vec/mode_3/mode_Ey_3.png " width="265"> <img src="./examples/modes_full_vec/mode_3/mode_Ez_3.png " width="265">

### Wavelength Sweep
<img src="./examples/modes_full_vec/wavelength_n_effs.png " width="400">

## Contributions
If you add functionality, please send me a pull request.

## Acknowledgments
This finite difference mode solver is based on a modified version of [EMpy](https://github.com/lbolla/EMpy).

Thank you to [Inna Krasnokutska](https://github.com/ikrasnokutska) for testing.
