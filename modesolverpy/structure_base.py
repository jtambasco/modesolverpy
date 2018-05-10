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
        interp_real = interpolate.interp2d(self.x, self.y, self.eps.real)
        interp_imag = interpolate.interp2d(self.x, self.y, self.eps.imag)
        interp = lambda x, y: interp_real(x, y) + 1.j*interp_imag(x, y)
        return interp

    @property
    def n_func(self):
        return interpolate.interp2d(self.x, self.y, self.n)

    def _add_triangular_sides(self, xy_mask, angle, y_top_right, y_bot_left,
                              x_top_right, x_bot_left, n_material):
        angle = np.radians(angle)
        trap_len = (y_top_right - y_bot_left) / np.tan(angle)
        num_x_iterations = round(trap_len/self.x_step)
        y_per_iteration = round(self.y_pts / num_x_iterations)

        lhs_x_start_index = int(x_bot_left/ self.x_step + 0.5)
        rhs_x_stop_index = int(x_top_right/ self.x_step + 1 + 0.5)

        for i, _ in enumerate(xy_mask):
            xy_mask[i][:lhs_x_start_index] = False
            xy_mask[i][lhs_x_start_index:rhs_x_stop_index] = True

            if i % y_per_iteration == 0:
                lhs_x_start_index -= 1
                rhs_x_stop_index += 1

        self.n[xy_mask] = n_material
        return self.n

    def add_material(self, x_bot_left, y_bot_left, x_top_right, y_top_right,
                     n_material, angle=0):

        x_mask = np.logical_and(x_bot_left<=self.x, self.x<=x_top_right)
        y_mask = np.logical_and(y_bot_left<=self.y, self.y<=y_top_right)

        xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
        self.n[xy_mask] = n_material

        if angle:
            self._add_triangular_sides(xy_mask, angle, y_top_right, y_bot_left,
                                       x_top_right, x_bot_left, n_material)

        return self.n

    def write_to_file(self, filename='material_index.dat', plot=True):
        path = os.path.dirname(sys.modules[__name__].__file__) + '/'

        with open(filename, 'w') as fs:
            for n_row in np.abs(self.n[::-1]):
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
        self._n = np.ones((self.y.size,self.x.size), 'complex_') * n_background

    @property
    def n(self):
        return self._n

class Slabs(_AbstractStructure):
    def __init__(self, wavelength, y_step, x_step, x_max, x_min=0.):
        _AbstractStructure.__init__(self)

        self._wl = wavelength
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

        if not callable(n_background):
            n_back = lambda wl: n_background
        else:
            n_back = n_background

        height_discretised = self.y_step*((height // self.y_step) + 1)

        y_min = self._next_start
        y_max = y_min + height_discretised
        self.slabs[name] = Slab(name, self.x_step, self.y_step, self.x_max,
                                y_max, self.x_min, y_min, n_back, self._wl)

        self.y_max = y_max
        self._next_start = y_min + height_discretised
        self.slab_count += 1

        return name

    def change_wavelength(self, wavelength):
        for name, slab in self.slabs.items():
            const_args = slab._const_args
            mat_args = slab._mat_params

            const_args[8] = wavelength

            s = Slab(*const_args)
            for mat_arg in mat_args:
                s.add_material(*mat_arg)

            self.slabs[name] = s

        self._wl = wavelength

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
                 n_background, wavelength):
        self._wl = wavelength
        self.name = name
        self.position = Slab.position
        Slab.position += 1

        Structure.__init__(self, x_step, y_step, x_max, y_max, x_min, y_min,
                           n_background(self._wl))

        self._const_args = [name, x_step, y_step, x_max, y_max, x_min, y_min, n_background, wavelength]
        self._mat_params = []

    def add_material(self, x_min, x_max, n, angle=0):
        self._mat_params.append([x_min, x_max, n, angle])

        if not callable(n):
            n_mat = lambda wl: n
        else:
            n_mat = n

        Structure.add_material(self, x_min, self.y_min, x_max, self.y_max, n_mat(self._wl), angle)
        return self.n

class StructureAni():
    def __init__(self, structure_xx, structure_yy, structure_zz,
                 structure_xy=None, structure_yx=None):
        self.xx = structure_yy
        self.yy = structure_xx
        self.zz = structure_zz

        if not structure_xy or not structure_yx:
            struct_dummy = Structure(self.xx.x_step, self.xx.y_step,
                                     self.xx.x_max, self.xx.y_max,
                                     self.xx.x_min, self.xx.y_min,
                                     n_background=0.)
            struct_dummy._wl = self.xx._wl

        if structure_xy:
            self.yx = structure_xy
        else:
            self.yx = struct_dummy

        if structure_yx:
            self.xy = structure_yx
        else:
            self.xy = struct_dummy

        assert self.xx._wl == self.xy._wl == self.yx._wl == \
               self.yy._wl == self.zz._wl

        self._wl = structure_xx._wl

        self.axes = (self.xx, self.xy, self.yx, self.yy, self.zz)
        self.axes_str = ('xx', 'xy', 'yx', 'yy', 'zz')

    @property
    def n(self):
        return [a.n for a in self.axes]

    @property
    def x_step(self):
        return self.xx.x_step

    @property
    def y_step(self):
        return self.xx.y_step

    @property
    def x_pts(self):
        return int((self.xx.x_max - self.xx.x_min) / self.xx.x_step + 1)

    @property
    def y_pts(self):
        return int((self.xx.y_max - self.xx.y_min) / self.xx.y_step)

    @property
    def x_ctr(self):
        return 0.5*(self.xx.x_max + self.xx.x_min)

    @property
    def y_ctr(self):
        return 0.5*(self.xx.y_max + self.xx.y_min)

    @property
    def xc(self):
        return 0.5*(self.xx.x[1:] + self.xx.x[:-1])

    @property
    def yc(self):
        return 0.5*(self.xx.y[1:] + self.xx.y[:-1])

    @property
    def xc_pts(self):
        return self.xx.x_pts - 1

    @property
    def yc_pts(self):
        return self.xx.y_pts - 1

    @property
    def xc_min(self):
        return self.xx.xc[0]

    @property
    def xc_max(self):
        return self.xx.xc[-1]

    @property
    def yc_min(self):
        return self.xx.yc[0]

    @property
    def yc_max(self):
        return self.xx.yc[-1]

    @property
    def x(self):
        if None not in (self.xx.x_min, self.xx.x_max, self.xx.x_step) and \
                self.xx.x_min != self.xx.x_max:
            x = np.arange(self.xx.x_min, self.xx.x_max+self.xx.x_step-self.xx.y_step*0.1, self.xx.x_step)
        else:
            x = np.array([])
        return x

    @property
    def y(self):
        if None not in (self.xx.y_min, self.xx.y_max, self.xx.y_step) and \
                self.xx.y_min != self.xx.y_max:
            y = np.arange(self.xx.y_min, self.xx.y_max-self.xx.y_step*0.1, self.xx.y_step)
        else:
            y = np.array([])
        return y

    @property
    def eps(self):
        eps_ani = [a.n**2 for a in self.axes]
        return eps_ani

    @property
    def eps_func(self):
        return lambda x,y: tuple(axis.eps_func(x,y) for axis in self.axes)

    @property
    def n_func(self):
        return lambda x,y: tuple(axis.n_func(x,y) for axis in self.axes)

    def write_to_file(self, filename='material_index.dat', plot=True):
        path = os.path.dirname(sys.modules[__name__].__file__) + '/'

        dir_plot = 'material_index/'
        if not os.path.exists(dir_plot):
            os.makedirs(dir_plot)

        for axis, name in zip(self.axes, self.axes_str):
            root, ext = os.path.splitext(filename)
            fn = dir_plot + root + '_'+ name + ext
            with open(fn, 'w') as fs:
                for n_row in np.abs(axis.n[::-1]):
                    n_str = ','.join([str(v) for v in n_row])
                    fs.write(n_str+'\n')

            if plot:
                filename_image_prefix, _ = os.path.splitext(fn)
                filename_image = filename_image_prefix + '.png'
                args = {
                    'title': 'Refractive Index Profile: %s' % name,
                    'x_pts': self.xx.x_pts,
                    'y_pts': self.xx.y_pts,
                    'x_min': self.xx.x_min,
                    'x_max': self.xx.x_max,
                    'y_min': self.xx.y_min,
                    'y_max': self.xx.y_max,
                    'filename_data': fn,
                    'filename_image': filename_image
                }
                gp.gnuplot(path+'structure.gpi', args, silent=False)

    def change_wavelength(self, wavelength):
        for axis in self.axes:
            if issubclass(type(axis), Slabs):
                axis.change_wavelength(wavelength)
        self.xx, self.xy, self.yx, self.yy, self.zz = self.axes
        self._wl = wavelength
