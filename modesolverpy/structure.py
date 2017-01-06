import copy
from .structure_base import *
import opticalmaterialspy as mat

class Updateable:
    def update(self, var, var_str, **kwargs):
        kwargs.update(self._vars)
        if var_str == 'step':
            kwargs['x_step'] = var
            kwargs['y_step'] = var
        else:
            kwargs[var_str] = var
        WgArray.__init__(**kwargs)

class RidgeWaveguide(Slabs, Updateable):
    def __init__(self, x_step, y_step, wg_height, wg_width, sub_height, sub_width,
                 clad_height, n_sub, n_wg, angle, trap_len = 0, n_clad=1., film_thickness='wg_height'):
        Slabs.__init__(self, y_step, x_step, sub_width)

        self.n_sub = n_sub
        self.n_clad = n_clad
        self.n_wg = n_wg

        self.add_slab(sub_height, n_sub)
        if film_thickness != 'wg_height' and film_thickness != wg_height:
            assert film_thickness > 0.
            self.add_slab(film_thickness-wg_height, n_wg)
        k = self.add_slab(wg_height, n_clad)

        if angle and not trap_len:
            trap_len = get_trap_length(angle, wg_height)
            trap_side = np.sqrt(trap_len**2 + wg_height**2)
            # print('trapside:', trap_side)
            # j = self.add_slab(trap_side, n_clad)
        else:
            trap_side = 0
        # self.slabs[k].add_material(self.x_ctr-wg_width-trap_len/2., self.x_ctr+wg_width+trap_len/2.,
        #                            n_wg)
        self.slabs[k].add_material(self.x_ctr-wg_width/2., self.x_ctr+wg_width/2.,
                                   n_wg, trap_side, trap_len)
        # self.slabs[j].add_material(self.x_ctr-wg_width-trap_len/2., self.x_ctr+wg_width+trap_len/2.,
        #                            n_wg)
        self.add_slab(clad_height, n_clad)

class WgArray(Slabs, Updateable):
    def __init__(self, x_step, y_step, wg_height, wg_widths, wg_gaps, sub_height,
                 sub_width, clad_height, n_sub, n_wg, angle, trap_len = 0, n_clad=1.):

        self._vars = locals()

        Slabs.__init__(self, y_step, x_step, sub_width)

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

        if angle:
            trap_len = get_trap_length(angle, wg_height)
            trap_side = np.sqrt(trap_len**2 + wg_height**2)
        # print(trap_len)
        k = self.add_slab(wg_height, n_clad)
        air_width_l_r = 0.5*(sub_width - np.sum(wg_widths) - np.sum(wg_gaps))
        position = air_width_l_r

        for wg_width, wg_gap in zip(wg_widths, wg_gaps_pad):
            self.slabs[k].add_material(position - trap_len, position+wg_width + trap_len, n_wg, trap_side, trap_len)
            print('position:',position)
            position += wg_width + wg_gap

        self.add_slab(clad_height, n_clad)

def get_trap_length(angle, wg_height):
    return wg_height / np.tan(angle)
