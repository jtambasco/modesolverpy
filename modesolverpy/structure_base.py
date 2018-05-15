import numpy as np
from scipy import interpolate
import gnuplotpy as gp
import os
import sys
import abc

class _AbstractStructure(metaclass=abc.ABCMeta):
    @abc.abstractproperty
    def n(self):
        '''
        np.array: A grid of refractive indices representing
        the refractive index profile of the structure.
        '''
        pass

    @property
    def x_pts(self):
        '''
        int: The number of grid points in x.
        '''
        return int((self.x_max - self.x_min) / self.x_step + 1)

    @property
    def y_pts(self):
        '''
        int: The number of grid points in y.
        '''
        return int((self.y_max - self.y_min) / self.y_step)

    @property
    def x_ctr(self):
        '''
        float: The centre distance in x.
        '''
        return 0.5*(self.x_max + self.x_min)

    @property
    def y_ctr(self):
        '''
        float: The centre distance in y
        '''
        return 0.5*(self.y_max + self.y_min)

    @property
    def xc(self):
        '''
        np.array: The centre points of the x points.
        '''
        return 0.5*(self.x[1:] + self.x[:-1])

    @property
    def yc(self):
        '''
        np.array: The centre points of the y points.
        '''
        return 0.5*(self.y[1:] + self.y[:-1])

    @property
    def xc_pts(self):
        '''
        int: The number of points in `xc`.
        '''
        return self.x_pts - 1

    @property
    def yc_pts(self):
        '''
        int: The number of points in `yc`.
        '''
        return self.y_pts - 1

    @property
    def xc_min(self):
        '''
        float: The minimum value of `xc`.
        '''
        return self.xc[0]

    @property
    def xc_max(self):
        '''
        float: The maximum value of `xc`.
        '''
        return self.xc[-1]

    @property
    def yc_min(self):
        '''
        float: The minimum value of `yc`.
        '''
        return self.yc[0]

    @property
    def yc_max(self):
        '''
        float: The maximum value of `yc`.
        '''
        return self.yc[-1]

    @property
    def x(self):
        '''
        np.array: The grid points in x.
        '''
        if None not in (self.x_min, self.x_max, self.x_step) and \
                self.x_min != self.x_max:
            x = np.arange(self.x_min, self.x_max+self.x_step-self.y_step*0.1, self.x_step)
        else:
            x = np.array([])
        return x

    @property
    def y(self):
        '''
        np.array: The grid points in y.
        '''
        if None not in (self.y_min, self.y_max, self.y_step) and \
                self.y_min != self.y_max:
            y = np.arange(self.y_min, self.y_max-self.y_step*0.1, self.y_step)
        else:
            y = np.array([])
        return y

    @property
    def eps(self):
        '''
        np.array: A grid of permittivies representing
        the permittivity profile of the structure.
        '''
        return self.n**2

    @property
    def eps_func(self):
        '''
        function: a function that when passed a `x` and `y` values,
            returns the permittivity profile of the structure,
            interpolating if necessary.
        '''
        interp_real = interpolate.interp2d(self.x, self.y, self.eps.real)
        interp_imag = interpolate.interp2d(self.x, self.y, self.eps.imag)
        interp = lambda x, y: interp_real(x, y) + 1.j*interp_imag(x, y)
        return interp

    @property
    def n_func(self):
        '''
        function: a function that when passed a `x` and `y` values,
            returns the refractive index profile of the structure,
            interpolating if necessary.
        '''
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

    def _add_material(self, x_bot_left, y_bot_left, x_top_right, y_top_right,
                     n_material, angle=0):
        '''
        A low-level function that allows writing a rectangle refractive
        index profile to a `Structure`.

        Args:
            x_bot_left (float): The bottom-left x-coordinate of the
                rectangle.
            y_bot_left (float): The bottom-left y-coordinate of the
                rectangle.
            x_top_right (float): The top-right x-coordinate of the
                rectangle.
            y_top_right (float): The top-right y-coordinate of the
                rectangle.
            n_material (float): The refractive index of the points
                encompassed by the defined rectangle.
            angle (float): The angle in degrees of the sidewalls
                of the defined rectangle.  Default is 0.  This
                is useful for creating a ridge with angled
                sidewalls.
        '''
        x_mask = np.logical_and(x_bot_left<=self.x, self.x<=x_top_right)
        y_mask = np.logical_and(y_bot_left<=self.y, self.y<=y_top_right)

        xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
        self.n[xy_mask] = n_material

        if angle:
            self._add_triangular_sides(xy_mask, angle, y_top_right, y_bot_left,
                                       x_top_right, x_bot_left, n_material)

        return self.n

    def write_to_file(self, filename='material_index.dat', plot=True):
        '''
        Write the refractive index profile to file.

        Args:
            filename (str): The nominal filename the refractive
                index data should be saved to.
            plot (bool): `True` if plots should be generates,
                otherwise `False`.  Default is `True`.
        '''
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
    '''
    Class to implement device refractive index
    profile cross-section designs.

    :class:`Slabs` is a collection of :class:`Slab` objects.  Each
    slab has a fixed height (usually less than the
    maximum height of the desired simulation window),
    and is as wide as the simulation window.

    :class:`Slabs` objects can be index using `[name]` to return
    the various :class:`Slab` objects.  The bottom slab is
    returned first and so on up to the top slab.

    .. image:: ../images/slabs.svg
        :width: 200%

    Args:
        wavelength (float): The wavelength the structure
            operates at.
        y_step (float): The step in y.
        x_step (float): The step in x.
        x_max (float): The maximum x-value.
        x_min (float): The minimum x-value. Default is 0.

    Attributes:
        slabs (dict): The key is the name of the slab,
            and the value is the :class:`Slab` object.
        slab_count (int): The number of :class:`Slab` objects
            added so far.
    '''
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
        '''
        Creates and adds a :class:`Slab` object.

        Args:
            height (float): Height of the slab.
            n_background (float): The nominal refractive
                index of the slab.  Default is 1 (air).

        Returns:
            str: The name of the slab.
        '''
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
        '''
        Changes the wavelength of the structure.

        This will affect the mode solver and potentially
        the refractive indices used (provided functions
        were provided as refractive indices).

        Args:
            wavelength (float): The new wavelength.
        '''
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
        '''
        np.array: The refractive index profile matrix
        of the current slab.
        '''
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
    '''
    A :class:`Slab` represents a horizontal slice of
    the refractive index profile.

    A :class:`Slabs` object composes many :class:`Slab` objects.
    The more :class:`Slab` are added, the more horizontal
    slices are added.  A :class:`Slab` has a chosen fixed
    height, and a background (nominal) refractive
    index.  A slab can then be customised to include
    a desired design.

    Args:
        name (str): The name of the slab.
        x_step (float): The step in x.
        y_step (float): The step in y.
        x_max (float): The maximum x-value.
        y_max (float): The maximum y-value.
        x_min (float): The minimum x-value.
        y_min (float): The minimum x-value.
        n_background (float): The nominal refractive
            index.
        wavelength (float): The wavelength the structure
            operates at.

    Attributes:
        name (str): The name of the :class:`Slab` object.
        position (int): A unique identifier for the
        :class:`Slab` object.
    '''
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
        '''
        Add a refractive index between two x-points.

        Args:
            x_min (float): The start x-point.
            x_max (float): The stop x-point.
            n (float, function):  Refractive index between
                 `x_min` and `x_max`.  Either a constant (`float`), or
                 a function that accepts one parameters, the
                 wavelength, and returns a float of the refractive
                 index.  This is useful when doing wavelength
                 sweeps and solving for the group velocity.  The
                 function provided could be a Sellmeier equation.
            angle (float): Angle in degrees of the slope of the
                sidewalls at `x_min` and `x_max`.  This is useful
                for defining a ridge with angled sidewalls.
        '''
        self._mat_params.append([x_min, x_max, n, angle])

        if not callable(n):
            n_mat = lambda wl: n
        else:
            n_mat = n

        Structure._add_material(self, x_min, self.y_min, x_max, self.y_max, n_mat(self._wl), angle)
        return self.n

class StructureAni():
    r"""
    Anisottropic structure object.

    This is used with the fully-vectorial simulation when
    an anisotropic material is being used.

    The form of the refractive index is

    .. math::

        n = \begin{bmatrix}
                n_{xx} & n_{xy} & 0 \\
                n_{yx} & n_{yy} & 0 \\
                0      & 0      & n_{zz}
            \end{bmatrix}.

    Args:
        structure_xx (Structure): The structure with refractive
            index, :math:`n_{xx}`.
        structure_yy (Structure): The structure with refractive
            index, :math:`n_{yy}`.  Presumably the same structure
            as `structure_xx`, but with different refractive index
            parameters.
        structure_zz (Structure): The structure with refractive
            index, :math:`n_{zz}`.  Presumably the same structure
            as `structure_xx`, but with different refractive index
            parameters.
        structure_xy (None, Structure): The structure with refractive
            index, :math:`n_{yx}`.  Presumably the same structure
            as `structure_xx`, but with different refractive index
            parameters.  Default is `None`.
        structure_yx (None, Structure): The structure with refractive
            index, :math:`n_{yx}`.  Presumably the same structure
            as `structure_xx`, but with different refractive index
            parameters.  Default is `None`.
    """
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
        '''
        Write the refractive index profile to file.

        Args:
            filename (str): The nominal filename the refractive
                index data should be saved to.
            plot (bool): `True` if plots should be generates,
                otherwise `False`.  Default is `True`.
        '''
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
        '''
        Changes the wavelength of the structure.

        This will affect the mode solver and potentially
        the refractive indices used (provided functions
        were provided as refractive indices).

        Args:
            wavelength (float): The new wavelength.
        '''
        for axis in self.axes:
            if issubclass(type(axis), Slabs):
                axis.change_wavelength(wavelength)
        self.xx, self.xy, self.yx, self.yy, self.zz = self.axes
        self._wl = wavelength
