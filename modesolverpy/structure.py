import copy
from . import structure_base as sb
from .structure_base import use_gnuplot
from .structure_base import use_matplotlib
import opticalmaterialspy as mat
import numpy as np

class RidgeWaveguide(sb.Slabs):
    '''
    A general ridge waveguide structure.

    Args:
        wavelength (float): Wavelength the structure should
            operate at.
        x_step (float): The grid step in x that the structure
            is created on.
        y_step (float): The grid step in y that the structure
            is created on.
        wg_height (float): The height of the ridge.
        wg_width (float): The width of the ridge.
        sub_height (float): The thickness of the substrate.
        sub_width (float): The width of the substrate.
        clad_height (float): The thickness of the cladding.
        n_sub (float, function): Refractive index of the
            substrate.  Either a constant (`float`), or
            a function that accepts one parameters, the
            wavelength, and returns a float of the refractive
            index.  This is useful when doing wavelength
            sweeps and solving for the group velocity.  The
            function provided could be a Sellmeier equation.
        n_wg (float, function): Refractive index of the
            waveguide.  Either a constant (`float`), or
            a function that accepts one parameters, the
            wavelength, and returns a float of the refractive
            index.  This is useful when doing wavelength
            sweeps and solving for the group velocity.  The
            function provided could be a Sellmeier equation.
        angle (float): The angle of the sidewall [degrees] of
            the waveguide.  Default is 0 degrees (vertical
            sidewalls).
        n_clad (float, function): Refractive index of the
            cladding.  Either a constant (`float`), or
            a function that accepts one parameters, the
            wavelength, and returns a float of the refractive
            index.  This is useful when doing wavelength
            sweeps and solving for the group velocity.  The
            function provided could be a Sellmeier equation.
            Default is air.
        film_thickness (float, str): The thickness of the
            film the waveguide is on.  If the waveguide
            is a true ridge (fully etched), then the film thickness
            can be set to 'wg_height', otherwise the waveguide
            is a rib waveguide, and a float should be given
            specifying the thickness of the film.

    '''
    def __init__(self, wavelength, x_step, y_step, wg_height, wg_width, sub_height, sub_width,
                 clad_height, n_sub, n_wg, angle=0, n_clad=mat.Air().n(),
                 film_thickness='wg_height'):
        sb.Slabs.__init__(self, wavelength, y_step, x_step, sub_width)

        self.n_sub = n_sub
        self.n_clad = n_clad
        self.n_wg = n_wg

        self.add_slab(sub_height, n_sub)
        if film_thickness != 'wg_height' and film_thickness != wg_height:
            assert film_thickness > 0., 'Film must have some thickness to it.'
            assert wg_height <= film_thickness, 'Waveguide can\'t be thicker than the film.'
            self.add_slab(film_thickness-wg_height, n_wg)
        k = self.add_slab(wg_height, n_clad)

        self.slabs[k].add_material(self.x_ctr-wg_width/2., self.x_ctr+wg_width/2.,
                                   n_wg, angle)

        self.add_slab(clad_height, n_clad)

class WgArray(sb.Slabs):
    def __init__(self, wavelength, x_step, y_step, wg_height, wg_widths, wg_gaps, sub_height,
                 sub_width, clad_height, n_sub, n_wg, angle=0, n_clad=mat.Air().n()):
        sb.Slabs.__init__(self, wavelength, y_step, x_step, sub_width)

        try:
            iter(wg_gaps)
        except TypeError:
            wg_gaps = [wg_gaps]

        try:
            assert len(wg_widths) == len(wg_gaps)+1
        except TypeError:
            wg_widths = [wg_widths for _ in wg_gaps]

        wg_gaps_pad = copy.copy(wg_gaps)
        wg_gaps_pad.append(0.)

        self.add_slab(sub_height, n_sub)

        k = self.add_slab(wg_height, n_clad)
        air_width_l_r = 0.5*(sub_width - np.sum(wg_widths) - np.sum(wg_gaps))
        position = air_width_l_r

        for wg_width, wg_gap in zip(wg_widths, wg_gaps_pad):
            self.slabs[k].add_material(position, position + wg_width, n_wg, angle)

            position += wg_width + wg_gap

        self.add_slab(clad_height, n_clad)
