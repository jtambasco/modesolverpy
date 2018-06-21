import abc
import os
import sys
import copy
import tqdm
import numpy as np
import gnuplotpy as gp
from . import _mode_solver_lib as ms
from . import _analyse as anal
from . import structure_base as stb

class _ModeSolver(metaclass=abc.ABCMeta):
    def __init__(self, n_eigs, tol=0., boundary='0000',
                 mode_profiles=True, initial_mode_guess=None,
                 n_eff_guess=None):
        self._n_eigs = int(n_eigs)
        self._tol = tol
        self._boundary = boundary
        self._mode_profiles = mode_profiles
        self._initial_mode_guess = initial_mode_guess
        self._n_eff_guess = n_eff_guess

        self.n_effs = None
        self.modes = None
        self.mode_types = None
        self.overlaps = None

        self._path = os.path.dirname(sys.modules[__name__].__file__) + '/'

    @abc.abstractproperty
    def _modes_directory(self):
        pass

    @abc.abstractmethod
    def _solve(self, structure, wavelength):
        pass

    def solve(self, structure):
        '''
        Find the modes of a given structure.

        Args:
            structure (Structure): The target structure to solve
                for modes.

        Returns:
            dict: The 'n_effs' key gives the effective indices
            of the modes.  The 'modes' key exists of mode
            profiles were solved for; in this case, it will
            return arrays of the mode profiles.
        '''
        return self._solve(structure, structure._wl)

    def solve_sweep_structure(self, structures, sweep_param_list,
                              filename='structure_n_effs.dat', plot=True):
        '''
        Find the modes of many structures.

        Args:
            structures (list): A list of `Structures` to find the modes
                of.
            sweep_param_list (list): A list of the parameter-sweep sweep
                that was used.  This is for plotting purposes only.
            filename (str): The nominal filename to use when saving the
                effective indices.  Defaults to 'structure_n_effs.dat'.
            plot (bool): `True` if plots should be generates,
                otherwise `False`.  Default is `True`.

        Returns:
            list: A list of the effective indices found for each structure.
        '''
        n_effs = []
        for s in tqdm.tqdm(structures, ncols=70):
            self.solve(s)
            n_effs.append(np.real(self.n_effs))

        if filename:
            self._write_n_effs_to_file(n_effs,
                                       self._modes_directory+filename,
                                       sweep_param_list)
            if plot:
                self._plot_n_effs(self._modes_directory+filename,
                                  'Structure number',
                                  'n_{effs} vs structure')

        return n_effs

    def solve_sweep_wavelength(self, structure, wavelengths, filename='wavelength_n_effs.dat',
                               plot=True):
        '''
        Solve for the effective indices of a fixed structure at
        different wavelengths.

        Args:
            structure (Slabs): The target structure to solve
                for modes.
            wavelengths (list): A list of wavelengths to sweep
                over.
            filename (str): The nominal filename to use when saving the
                effective indices.  Defaults to 'wavelength_n_effs.dat'.
            plot (bool): `True` if plots should be generates,
                otherwise `False`.  Default is `True`.

        Returns:
            list: A list of the effective indices found for each wavelength.
        '''
        n_effs = []
        for w in tqdm.tqdm(wavelengths, ncols=70):
            structure.change_wavelength(w)
            self.solve(structure)
            n_effs.append(np.real(self.n_effs))

        if filename:
            self._write_n_effs_to_file(n_effs,
                                       self._modes_directory+filename,
                                       wavelengths)
            if plot:
                self._plot_n_effs(self._modes_directory+filename,
                                  'Wavelength',
                                  'n_{effs} vs wavelength')

        return n_effs

    def solve_ng(self, structure, wavelength_step=0.01, filename='ng.dat'):
        r'''
        Solve for the group index, :math:`n_g`, of a structure at a particular
        wavelength.

        Args:
            structure (Structure): The target structure to solve
                for modes.
            wavelength_step (float): The step to take below and
                above the nominal wavelength.  This is used for
                approximating the gradient of :math:`n_\mathrm{eff}`
                at the nominal wavelength.  Default is 0.01.
            filename (str): The nominal filename to use when saving the
                effective indices.  Defaults to 'wavelength_n_effs.dat'.

        Returns:
            list: A list of the group indices found for each mode.
        '''
        wl_nom = structure._wl

        self.solve(structure)
        n_ctrs = self.n_effs

        structure.change_wavelength(wl_nom-wavelength_step)
        self.solve(structure)
        n_bcks = self.n_effs

        structure.change_wavelength(wl_nom+wavelength_step)
        self.solve(structure)
        n_frws = self.n_effs

        n_gs = []
        for n_ctr, n_bck, n_frw in zip(n_ctrs, n_bcks, n_frws):
            n_gs.append(n_ctr - wl_nom*(n_frw-n_bck)/(2*wavelength_step))

        if filename:
            with open(self._modes_directory+filename, 'w') as fs:
                fs.write('# Mode idx, Group index\n')
                for idx, n_g in enumerate(n_gs):
                    fs.write('%i,%.3f\n' % (idx, np.round(n_g.real,3)))

        return n_gs

    def _get_mode_filename(self, field_name, mode_number, filename):
        filename_prefix, filename_ext = os.path.splitext(filename)
        filename_mode = filename_prefix + '_' + field_name + \
            '_' + str(mode_number) + filename_ext
        return filename_mode

    def _write_n_effs_to_file(self, n_effs, filename, x_vals=None):
        with open(filename, 'w') as fs:
            for i, n_eff in enumerate(n_effs):
                if x_vals is not None:
                    line_start = str(x_vals[i])+','
                else:
                    line_start = ''
                line = ','.join([str(np.round(n,3)) for n in n_eff])
                fs.write(line_start+line+'\n')
        return n_effs

    def _write_mode_to_file(self, mode, filename):
        with open(filename, 'w') as fs:
            for e in mode[::-1]:
                e_str = ','.join([str(v) for v in e])
                fs.write(e_str+'\n')
        return mode

    def _plot_n_effs(self, filename_n_effs, xlabel, title):
        args = {
            'titl': title,
            'xlab': xlabel,
            'ylab': 'n_{eff}',
            'filename_data': filename_n_effs,
            'filename_image': None,
            'num_modes': len(self.modes)
        }

        filename_image_prefix, _ = os.path.splitext(filename_n_effs)
        filename_image = filename_image_prefix + '.png'
        args['filename_image'] = filename_image

        gp.gnuplot(self._path+'n_effs.gpi', args, silent=False)
        gp.trim_pad_image(filename_image)

        return args

    def _plot_mode(self, field_name, mode_number, filename_mode, n_eff=None,
                   subtitle='', e2_x=0., e2_y=0., ctr_x=0., ctr_y=0.,
                   area=None, wavelength=None):
        fn = field_name[0] + '_{' + field_name[1:] + '}'
        title = 'Mode %i |%s| Profile' % (mode_number, fn)
        if n_eff:
            title += ', n_{eff}: ' + '{:.3f}'.format(n_eff.real)
        if wavelength:
            title += ', λ = %s ' % '{:.3f} µm'.format(wavelength)
        if area:
            title += ', A_%s: ' % field_name[1] + '{:.1f}\%'.format(area)

        if subtitle:
            title += '\n{/*0.7 %s}' % subtitle

        args = {
            'title': title,
            'x_pts': self._structure.xc_pts,
            'y_pts': self._structure.yc_pts,
            'x_min': self._structure.xc_min,
            'x_max': self._structure.xc_max,
            'y_min': self._structure.yc_min,
            'y_max': self._structure.yc_max,
            'x_step': self._structure.x_step,
            'y_step': self._structure.y_step,
            'filename_data': filename_mode,
            'filename_image': None,
            'e2_x': e2_x,
            'e2_y': e2_y,
            'ctr_x': ctr_x,
            'ctr_y': ctr_y
        }

        filename_image_prefix, _ = os.path.splitext(filename_mode)
        filename_image = filename_image_prefix + '.png'
        args['filename_image'] = filename_image

        gp.gnuplot(self._path+'mode.gpi', args)
        gp.trim_pad_image(filename_image)

        return args

class ModeSolverSemiVectorial(_ModeSolver):
    '''
    A semi-vectorial mode solver object used to
    setup and run a mode solving simulation.

    Args:
        n_eigs (int): The number of eigen-values to solve for.
        tol (float): The precision of the eigen-value/eigen-vector
            solver.  Default is 0.001.
        boundary (str): The boundary conditions to use.
            This is a string that identifies the type of boundary conditions applied.
            The following options are available: 'A' - Hx is antisymmetric, Hy is symmetric,
            'S' - Hx is symmetric and, Hy is antisymmetric, and '0' - Hx and Hy are zero
            immediately outside of the boundary.
            The string identifies all four boundary conditions, in the order:
            North, south, east, west. For example, boundary='000A'. Default is '0000'.
        mode_profiles (bool): `True if the the mode-profiles should be found, `False`
            if only the effective indices should be found.
        initial_mode_guess (list): An initial mode guess for the modesolver.
        semi_vectorial_method (str): Either 'Ex' or 'Ey'.  If 'Ex', the mode solver
            will only find TE modes (horizontally polarised to the simulation window),
            if 'Ey', the mode solver will find TM modes (vertically polarised to the
            simulation window).
    '''
    def __init__(self, n_eigs, tol=0.001, boundary='0000',
                 mode_profiles=True, initial_mode_guess=None,
                 semi_vectorial_method='Ex'):
        self._semi_vectorial_method = semi_vectorial_method
        _ModeSolver.__init__(self, n_eigs, tol, boundary,
                             mode_profiles, initial_mode_guess)

    @property
    def _modes_directory(self):
        modes_directory = './modes_semi_vec/'
        if not os.path.exists(modes_directory):
            os.mkdir(modes_directory)
        _modes_directory = modes_directory
        return _modes_directory

    def _solve(self, structure, wavelength):
        self._structure = structure
        self._ms = ms._ModeSolverSemiVectorial(wavelength,
                                               structure,
                                               self._boundary,
                                               self._semi_vectorial_method)
        self._ms.solve(self._n_eigs, self._tol, self._mode_profiles,
                       initial_mode_guess=self._initial_mode_guess)
        self.n_effs = self._ms.neff

        r = {'n_effs': self.n_effs}

        if self._mode_profiles:
            r['modes'] = self._ms.modes
            self._ms.modes[0] = np.real(self._ms.modes[0])
            self._initial_mode_guess = np.real(self._ms.modes[0])

        self.modes = self._ms.modes

        return r

    def write_modes_to_file(self, filename='mode.dat', plot=True, analyse=True):
        '''
        Writes the mode fields to a file and optionally plots them.

        Args:
            filename (str): The nominal filename to use for the saved
                data.  The suffix will be automatically be changed to
                identifiy each mode number.  Default is 'mode.dat'
            plot (bool): `True` if plots should be generates,
                otherwise `False`.  Default is `True`.
            analyse (bool): `True` if an analysis on the fundamental
                mode should be performed.  The analysis adds to the
                plot of the fundamental mode the power mode-field
                diameter (MFD) and marks it on the output, and it
                marks with a cross the maximum E-field value.
                Default is `True`.

        Returns:
            dict: A dictionary containing the effective indices
            and mode field profiles (if solved for).
        '''
        modes_directory = './modes_semi_vec/'
        if not os.path.isdir(modes_directory):
            os.mkdir(modes_directory)
        filename = modes_directory + filename

        for i, mode in enumerate(self._ms.modes):
            filename_mode = self._get_mode_filename(self._semi_vectorial_method,
                                                    i, filename)
            self._write_mode_to_file(np.real(mode), filename_mode)

            if plot:
                if i == 0 and analyse:
                    A, centre, sigma_2 = anal.fit_gaussian(self._structure.xc,
                                                           self._structure.yc,
                                                           np.abs(mode))
                    subtitle = ('E_{max} = %.3f, (x_{max}, y_{max}) = (%.3f, %.3f), MFD_{x} = %.3f, '
                                'MFD_{y} = %.3f') % (A, centre[0], centre[1], sigma_2[0], sigma_2[1])
                    self._plot_mode(self._semi_vectorial_method, i, filename_mode,
                                    self.n_effs[i], subtitle, sigma_2[0], sigma_2[1],
                                    centre[0], centre[1], wavelength=self._structure._wl)
                else:
                    self._plot_mode(self._semi_vectorial_method, i, filename_mode,
                                    self.n_effs[i], wavelength=self._structure._wl)

        return self.modes

class ModeSolverFullyVectorial(_ModeSolver):
    '''
    A fully-vectorial mode solver object used to
    setup and run a mode solving simulation.

    Args:
        n_eigs (int): The number of eigen-values to solve for.
        tol (float): The precision of the eigen-value/eigen-vector
            solver.  Default is 0.001.
        boundary (str): The boundary conditions to use.
            This is a string that identifies the type of boundary conditions applied.
            The following options are available: 'A' - Hx is antisymmetric, Hy is symmetric,
            'S' - Hx is symmetric and, Hy is antisymmetric, and '0' - Hx and Hy are zero
            immediately outside of the boundary.
            The string identifies all four boundary conditions, in the order:
            North, south, east, west. For example, boundary='000A'. Default is '0000'.
        initial_mode_guess (list): An initial mode guess for the modesolver.
        initial_n_eff_guess (list): An initial effective index guess for the modesolver.
    '''
    def __init__(self, n_eigs, tol=0.001, boundary='0000',
                 initial_mode_guess=None, n_eff_guess=None):
        self.n_effs_te = None
        self.n_effs_tm = None
        _ModeSolver.__init__(self, n_eigs, tol, boundary,
                             False, initial_mode_guess,
                             n_eff_guess)

    @property
    def _modes_directory(self):
        modes_directory = './modes_full_vec/'
        if not os.path.exists(modes_directory):
            os.mkdir(modes_directory)
        _modes_directory = modes_directory
        return _modes_directory

    def _solve(self, structure, wavelength):
        self._structure = structure
        self._ms = ms._ModeSolverVectorial(wavelength,
                                           structure,
                                           self._boundary)
        self._ms.solve(self._n_eigs, self._tol, self._n_eff_guess,
                       initial_mode_guess=self._initial_mode_guess)
        self.n_effs = self._ms.neff

        r = {'n_effs': self.n_effs}
        r['modes'] = self.modes = self._ms.modes

        self.overlaps = self._get_overlaps(self.modes)
        self.mode_types = self._get_mode_types()

        self._initial_mode_guess = None

        self.n_effs_te, self.n_effs_tm = self._sort_neffs(self._ms.neff)

        return r

    def _get_mode_types(self):
        mode_types = []
        labels = {0: 'qTE', 1: 'qTM', 2: 'qTE/qTM'}
        for overlap in self.overlaps:
            idx = np.argmax(overlap[0:3])
            mode_types.append((labels[idx], np.round(overlap[idx],2)))
        return mode_types

    def _sort_neffs(self, n_effs):
        mode_types = self._get_mode_types()

        n_effs_te = []
        n_effs_tm = []

        for mt, n_eff in zip(mode_types, n_effs):
            if mt[0] == 'qTE':
                n_effs_te.append(n_eff)
            elif mt[0] == 'qTM':
                n_effs_tm.append(n_eff)

        return n_effs_te, n_effs_tm

    def _get_overlaps(self, fields):
        mode_areas = []
        for mode in self._ms.modes:
            e_fields = (mode.fields['Ex'], mode.fields['Ey'], mode.fields['Ez'])
            h_fields = (mode.fields['Hx'], mode.fields['Hy'], mode.fields['Hz'])

            areas_e = [np.sum(np.abs(e)**2) for e in e_fields]
            areas_e /= np.sum(areas_e)
            areas_e *= 100

            areas_h = [np.sum(np.abs(h)**2) for h in h_fields]
            areas_h /= np.sum(areas_h)
            areas_h *= 100

            areas = areas_e.tolist()
            areas.extend(areas_h)
            mode_areas.append(areas)

        return mode_areas

    def write_modes_to_file(self, filename='mode.dat', plot=True,
                            fields_to_write=('Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz')):
        '''
        Writes the mode fields to a file and optionally plots them.

        Args:
            filename (str): The nominal filename to use for the saved
                data.  The suffix will be automatically be changed to
                identifiy each field and mode number.  Default is
                'mode.dat'
            plot (bool): `True` if plots should be generates,
                otherwise `False`.  Default is `True`.
            fields_to_write (tuple): A tuple of strings where the
                strings can be 'Ex', 'Ey', 'Ez', 'Hx', 'Hy' and 'Hz'
                defining what part of the mode should be saved and
                plotted.  By default, all six components are written
                and plotted.

        Returns:
            dict: A dictionary containing the effective indices
            and mode field profiles (if solved for).
        '''
        modes_directory = self._modes_directory

        # Mode info file.
        with open(modes_directory+'mode_info', 'w') as fs:
            fs.write('# Mode idx, Mode type, % in major direction, n_eff\n')
            for i, (n_eff, (mode_type, percentage)) in enumerate(zip(self.n_effs, self.mode_types)):
                mode_idx = str(i)
                line = '%s,%s,%.2f,%.3f' % (mode_idx, mode_type, percentage, n_eff.real)
                fs.write(line+'\n')

        # Mode field plots.
        for i, (mode, areas) in enumerate(zip(self._ms.modes, self.overlaps)):
            mode_directory = '%smode_%i/' % (modes_directory, i)
            if not os.path.isdir(mode_directory):
                os.mkdir(mode_directory)
            filename_full = mode_directory + filename

            for (field_name, field_profile), area in zip(mode.fields.items(), areas):
                if field_name in fields_to_write:
                    filename_mode = self._get_mode_filename(field_name, i, filename_full)
                    self._write_mode_to_file(np.real(field_profile),
                                             filename_mode)
                    if plot:
                        self._plot_mode(field_name, i, filename_mode, self.n_effs[i],
                                        area=area, wavelength=self._structure._wl)

        return self.modes
