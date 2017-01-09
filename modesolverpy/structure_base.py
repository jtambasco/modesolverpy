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


    def _nxt_piece(self, x_mask, y_mask, n_material, m, num_y = 0):

        x_done = False
        count = 0
        for i in range(0, len(x_mask)):
            if x_mask[i]:
                j = i
                count += 1

        if count < 2:
            x_mask[j] = False
            x_done = True

        for i in range(0,len(x_mask)):
            if i < len(x_mask) + 1 and x_mask[m * i] and x_mask[m * (i + 1)]:
                if x_mask[m * (i + 1)] and x_mask[m * (i + 1)]:
                    x_mask[m * i] = False
                    break

        i = 1
        g = 0
        found = False
        while i < len(y_mask) and not found:
            if i < len(y_mask) - 1 and y_mask[-i] and not y_mask[-(i + 1)]:
                while g + i < len(y_mask) and g < num_y + 1:
                    y_mask[-(i + g)] = True
                    found = True
                    g += 1
            i += 1

        if count > 1:
            xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
            self.n[xy_mask] = n_material
            self._nxt_piece(x_mask, y_mask, n_material, m, num_y)
        else:
            xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
            self.n[xy_mask] = n_material
            return self.n

    def _right_diagonal(self, x_top_right, x_mask, y_mask, trap_len, n_material):

        x_done = False
        iterations = trap_len/self.x_step

        i = 0
        g = 0
        while i < len(x_mask) and not x_done:
            if i * self.x_step > x_top_right + trap_len:
                x_done = True

            if i * self.x_step > x_top_right and i * self.x_step < x_top_right + trap_len and not x_done:
                for g in range(0, int(iterations) + 1):
                    x_mask[i] = True
                    i += 1
                    x_done = True
            i += 1

        num_x = 0
        for i in range(0, len(x_mask)):
            if x_mask[i]:
                num_x += 1

        num_y = len(y_mask)
        num_y_per_round = num_y / num_x
        num_y_per_round = int(num_y_per_round) + 1

        y_mask[-1] = True

        xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
        self.n[xy_mask] = n_material

        self._nxt_piece(x_mask, y_mask, n_material, -1, num_y_per_round)

        return self.n

    def _left_diagonal(self, x_bot_left, x_mask, y_mask, trap_len, n_material):

        x_done = False
        iterations = trap_len/self.x_step

        i = 0
        g = 0
        while i < len(x_mask) and not x_done:
            if i * self.x_step > x_bot_left:
                x_done = True

            if i * self.x_step > x_bot_left - trap_len and i*self.x_step < x_bot_left and not x_done:
                for g in range(0, int(iterations) + 1):
                    x_mask[i] = True
                    i += 1
                    x_done = True
            i += 1

        num_x = 0
        for i in range(0, len(x_mask)):
            if x_mask[i]:
                num_x += 1

        num_y = len(y_mask)
        num_y_per_round = num_y / num_x
        num_y_per_round = int(num_y_per_round) + 1

        y_mask[-1] = True

        xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
        self.n[xy_mask] = n_material

        self._nxt_piece(x_mask, y_mask, n_material, 1, num_y_per_round)
        return self.n

    def add_material(self, x_bot_left, y_bot_left, x_top_right, y_top_right,
                     n_material, trap_len = 0):

        x_mask = np.logical_and(x_bot_left<=self.x, self.x<=x_top_right)
        y_mask = np.logical_and(y_bot_left<=self.y, self.y<=y_top_right)

        xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
        self.n[xy_mask] = n_material

        if trap_len:
            for i in range(0, len(x_mask)):
                x_mask[i] = False
            for i in range(0, len(y_mask)):
                y_mask[i] = False

            self._left_diagonal(x_bot_left, x_mask, y_mask, trap_len, n_material)

            for i in range(0, len(x_mask)):
                x_mask[i] = False
            for i in range(0, len(y_mask)):
                y_mask[i] = False

            self._right_diagonal(x_top_right, x_mask, y_mask, trap_len, n_material)

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
