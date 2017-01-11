import numpy as np
from scipy import interpolate
import gnuplotpy as gp
import os
import sys
import abc

class _AbstractStructure(metaclass=abc.ABCMeta):
    @abc.abstractproperty
    def n(self):
        pass

    @property
    def x_pts(self):
        return int((self.x_max - self.x_min) / self.x_step + 1)

    @property
    def y_pts(self):
        return int((self.y_max - self.y_min) / self.y_step)

    @property
    def x_ctr(self):
        return 0.5*(self.x_max + self.x_min)

    @property
    def y_ctr(self):
        return 0.5*(self.y_max + self.y_min)

    @property
    def xc(self):
        return 0.5*(self.x[1:] + self.x[:-1])

    @property
    def yc(self):
        return 0.5*(self.y[1:] + self.y[:-1])

    @property
    def xc_pts(self):
        return self.x_pts - 1

    @property
    def yc_pts(self):
        return self.y_pts - 1

    @property
    def xc_min(self):
        return self.xc[0]

    @property
    def xc_max(self):
        return self.xc[-1]

    @property
    def yc_min(self):
        return self.yc[0]

    @property
    def yc_max(self):
        return self.yc[-1]

    @property
    def x(self):
        if None not in (self.x_min, self.x_max, self.x_step) and \
                self.x_min != self.x_max:
            x = np.arange(self.x_min, self.x_max+self.x_step-self.y_step*0.1, self.x_step)
        else:
            x = np.array([])
        return x

    @property
    def y(self):
        if None not in (self.y_min, self.y_max, self.y_step) and \
                self.y_min != self.y_max:
            y = np.arange(self.y_min, self.y_max-self.y_step*0.1, self.y_step)
        else:
            y = np.array([])
        return y

    @property
    def eps(self):
        return self.n**2

    @property
    def eps_func(self):
        return interpolate.interp2d(self.x, self.y, self.eps)

    @property
    def n_func(self):
        return interpolate.interp2d(self.x, self.y, self.n)

    def _add_triangular_sides(self, xy_mask, trap_len, y_top_right, y_bot_left, x_top_right, x_bot_left, n_material):

        num_x_iterations = round(trap_len/self.x_step)
        y_per_iteration = round(self.y_pts / num_x_iterations)

        lhs_x_start_index = round(x_bot_left/ self.x_step)
        rhs_x_stop_index = round(x_top_right/ self.x_step) + 1

        for i, _ in enumerate(xy_mask):
            xy_mask[i][0:lhs_x_start_index] = False
            xy_mask[i][lhs_x_start_index:rhs_x_stop_index] = True

            if i % y_per_iteration == 0:
                lhs_x_start_index -= 1
                rhs_x_stop_index += 1

        self.n[xy_mask] = n_material
        return self.n

    def add_material(self, x_bot_left, y_bot_left, x_top_right, y_top_right,
                     n_material, trap_len = 0):

        x_mask = np.logical_and(x_bot_left<=self.x, self.x<=x_top_right)
        y_mask = np.logical_and(y_bot_left<=self.y, self.y<=y_top_right)

        xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
        self.n[xy_mask] = n_material

        if trap_len:
            self._add_triangular_sides(xy_mask,trap_len, y_top_right, y_bot_left, x_top_right, x_bot_left, n_material)

        return self.n

    def write_to_file(self, filename='material_index.dat', plot=True):
        path = os.path.dirname(sys.modules[__name__].__file__) + '/'

        with open(filename, 'w') as fs:
            for n_row in self.n[::-1]:
                n_str = ','.join([str(v) for v in n_row])
                fs.write(n_str+'\n')

        if plot:
            filename_image_prefix, _ = os.path.splitext(filename)
            filename_image = filename_image_prefix + '.png'
            args = {
                'title': 'Refractive Index Profile',
                'x_pts': self.x_pts,
                'y_pts': self.y_pts,
                'x_min': self.x_min,
                'x_max': self.x_max,
                'y_min': self.y_min,
                'y_max': self.y_max,
                'filename_data': filename,
                'filename_image': filename_image
            }
            gp.gnuplot(path+'structure.gpi', args)

    def __str__(self):
        return self.n.__str__()

class Structure(_AbstractStructure):
    def __init__(self, x_step, y_step, x_max, y_max, x_min=0., y_min=0.,
                 n_background=1.):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.x_step = x_step
        self.y_step = y_step
        self.n_background = n_background
        self._n = np.ones((self.y.size,self.x.size)) * n_background

    @property
    def n(self):
        return self._n

class Slabs(_AbstractStructure):
    def __init__(self, y_step, x_step, x_max, x_min=0.):
        _AbstractStructure.__init__(self)

        self.x_min = x_min
        self.x_max = x_max
        self.x_step = x_step
        self.y_step = y_step
        self.y_min = 0

        self.slabs = {}
        self.slab_count = 0
        self._next_start = 0.

    def add_slab(self, height, n_background=1.):
        name = str(self.slab_count)

        height_discretised = self.y_step*((height // self.y_step) + 1)

        y_min = self._next_start
        y_max = y_min + height_discretised
        self.slabs[name] = Slab(name, self.x_step, self.y_step, self.x_max,
                                y_max, self.x_min, y_min, n_background)

        self.y_max = y_max
        self._next_start = y_min + height_discretised
        self.slab_count += 1

        return name

    @property
    def n(self):
        try:
            n_mat = self.slabs['0'].n
            for s in range(1, self.slab_count):
                n_mat = np.vstack((self.slabs[str(s)].n, n_mat))
        except KeyError:
            n_mat = None
        return n_mat

    def __getitem__(self, slab_name):
        return self.slabs[str(slab_name)]

class Slab(Structure):
    position = 0

    def __init__(self, name, x_step, y_step, x_max, y_max, x_min, y_min,
                 n_background):
        Structure.__init__(self, x_step, y_step, x_max, y_max, x_min, y_min,
                           n_background)
        self.name = name
        self.position = Slab.position
        Slab.position += 1

    def add_material(self, x_min, x_max, n, trap_len = 0):
        Structure.add_material(self, x_min, self.y_min, x_max, self.y_max, n, trap_len)
        return self.n

