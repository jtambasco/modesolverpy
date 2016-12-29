# modesolverpy
Photonic mode solver with a simple interface.

This finite difference mode solver is essentially a butchered version of [EMpy](https://github.com/lbolla/EMpy), with a simply interface to create typical waveguide structures.

This was written to simplify the simulation of structures that I fabricate as part of my research, and I thought others may benefit from it too.

## Installation
It is recommend to install `modesolverpy` either via:

### Pip:
```bash
pip3 install git+https://github.com/jtambasco/gnuplotpy.git
pip3 install git+https://github.com/jtambasco/opticalmaterialspy.git # optional
pip3 install git+https://github.com/jtambasco/modesolverpy.git 
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
* [optional] [opticalmaterialspy](https://github.com/jtambasco/opticalmaterialspy).

#### Other

* [Gnuplot](http://www.gnuplot.info/).

## Features
The main reasons to consider this library include:

* semi-vectorial and fully vectorial options,
* simple structure drawing,
* automated data saving and plotting via Gnuplot,
* some limited (at this stage) data processing (finding MFD of fundamental mode),
* easily extensible library, and
* works well with [opticalmaterialspy](https://github.com/jtambasco/opticalmaterialspy).

## Example

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

mode_solver = ms.ModeSolverSemiVectorial(2)
mode_solver.solve(structure, wavelength)
mode_solver.write_modes_to_file('example_modes_1.dat')
```

#### Structure
<img src="./examples/example_structure_1.png " width="400">

### Modes
<img src="./examples/example_modes_1_Ex_0.png " width="400"> <img src="./examples/example_modes_1_Ex_1.png " width="400">

## Contributions
If you add functionality, I'd be interested and would appreciate if you send me a pull request.
