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

        modes_directory = './modes_full_vec/'
        if not os.path.exists(modes_directory):
            os.mkdir(modes_directory)
        self._modes_directory = modes_directory

    @abc.abstractmethod
    def _solve(self, structure, wavelength):
        pass

    def solve(self, structure):
        return self._solve(structure, structure._wl)

    def solve_sweep_structure(self, structures, sweep_param_list,
                              filename='structure_n_effs.dat', plot=True):
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

    def solve_ng(self, structure, wavelength, wavelength_step=0.01, filename='ng.dat'):
        self.solve(structure)
        n_ctrs = self.n_effs

        structure.change_wavelength(wavelength-wavelength_step)
        self.solve(structure)
        n_bcks = self.n_effs

        structure.change_wavelength(wavelength+wavelength_step)
        self.solve(structure)
        n_frws = self.n_effs

        n_gs = []
        for n_ctr, n_bck, n_frw in zip(n_ctrs, n_bcks, n_frws):
            n_gs.append(n_ctr - wavelength*(n_frw-n_bck)/(2*wavelength_step))

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
    def __init__(self, n_eigs, tol=0., boundary='0000',
                 mode_profiles=True, initial_mode_guess=None,
                 semi_vectorial_method='Ex'):
        self._semi_vectorial_method = semi_vectorial_method
        _ModeSolver.__init__(self, n_eigs, tol, boundary,
                             mode_profiles, initial_mode_guess)

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
    def __init__(self, n_eigs, tol=0.001, boundary='0000',
                 initial_mode_guess=None, n_eff_guess=None):
        _ModeSolver.__init__(self, n_eigs, tol, boundary,
                             False, initial_mode_guess,
                             n_eff_guess)

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

        return r

    def _get_mode_types(self):
        mode_types = []
        labels = {0: 'qTE', 1: 'qTM', 2: 'qTE/qTM'}
        for overlap in self.overlaps:
            idx = np.argmax(overlap[0:3])
            mode_types.append((labels[idx], np.round(overlap[idx],2)))
        return mode_types

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
