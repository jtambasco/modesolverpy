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


    def _nxt_piece(self, x_mask, y_mask, n_material, x_mask_ref, y_mask_ref, left):

        found = False
        x_done = False
        y_done = False

        if left:
            i = 0
            while i < len(x_mask) and not found:
                if x_mask[i]:
                    found = True
                    x_mask[i] = False
                i += 1

            count = 0
            for this in x_mask:
                if this == True:
                    count += 1

            i = 1
            done = False
            while i < len(y_mask) and not done:
                while i < len(y_mask) and y_mask[-i] and y_mask_ref[-i]:
                    i += 1
                while i < len(y_mask) and not y_mask_ref[-i]:
                    y_mask_ref[-i] = True
                    y_mask[-i] = True
                    i += 1
                    done = True
                if y_mask_ref[-i]:
                    y_mask[-i] = True

                i += 1

            i = 1
            while i < len(y_mask):
                if y_mask_ref[i] == y_mask[i]:
                    y_done = True
                else:
                    y_done = False
                    break
                i += 1
            # print('x_mask', x_mask)
            # print('y_mask', y_mask)

            if not y_done or count > 1:
                xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
                self.n[xy_mask] = n_material
                self._nxt_piece(x_mask, y_mask, n_material, x_mask_ref, y_mask_ref, left)
            else:
                xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
                self.n[xy_mask] = n_material
                return self.n

    def _right_diagonal(self, x_bot_left, y_bot_left, x_top_right, y_top_right,
                     n_material, trap_side = 0, trap_len = 0):

        r_x_max = []
        r_x_min = x_top_right

        iterations = trap_len/self.x_step
        discrete_step = trap_len / iterations

        trap_range = y_top_right - y_bot_left

        y_max = []
        y_max.append(self.y_max)
        diff = iterations - int(iterations)


        # print('r_x_max', r_x_max)

        new_y = []

        done = False
        i = 1

        print("\nBefore")
        print ('y_max', y_max)
        print ('self.y', self.y)

        while not len(y_max) == len(self.y):
            val = y_max[-1] - (trap_range / iterations)

            while i < len(self.y) and val - self.y[i] > self.x_step and not len(self.y) == len(y_max) + 1:
                # print('comparing:', val, 'against ', self.y[i])
                y_max.append(self.y[i] - 0.0000001)
                i += 1
            if val > y_top_right:
                y_max.append(y_top_right)
            else:
                y_max.append(val)
                new_y.append(i)
                i+=1

        print("\nAfter")
        print ('y_max', y_max)
        print ('self.y', self.y)


        # print('y_max', y_max)

        r_x_max = self.x + 1
        ref_x = self.x
        done = False
        i = 0
        # print("\nBEFORE:")
        # print('ref_x',ref_x)
        # print('\nr_x_max', r_x_max)

        while i < len(r_x_max):
            if r_x_max[i] - 1 < x_top_right + trap_len and r_x_max[i] - 1 > x_top_right and not done:
                for g in range(0, int(iterations)):
                    # print("BURRITO")
                    r_x_max[i] = (x_top_right + trap_len) - g * self.x_step
                    i += 1
                    done = True

                if done:
                    r_x_max[i] = x_top_right + trap_len
                    # if l_x_max % self.x_step != 0:
                    #     l_x_max += 0.01
                    ref_x[i] = x_top_right + trap_len - 0.0001
                    break
            i += 1

        # print('i:', i)
        # i += 1
        while done and i < len(r_x_max):
            r_x_max[i] -= 1.1
            i += 1

        print("\nAFTER:")
        print('ref_x',ref_x)
        print('\nr_x_max', r_x_max)


        done = False
        i = 0
        while i < len(y_max) - 1:
            if y_max[i] > self.y[i]:
                while g < int(iterations) and not done:
                    while y_max[i + 1] < self.y[i + 1]:
                        y_max[i + 1] = self.y[i + 1] + 0.01
                        i += 1
                        done = True
            i += 1

        y_top_right = y_max

        x_mask_ref = []
        y_mask_ref = []

        x_mask = np.logical_and(r_x_min<=ref_x, ref_x<=r_x_max)
        y_mask = np.logical_and(y_bot_left<=self.y, self.y<=y_top_right)

        for i in range(0, len(x_mask)):
            x_mask_ref.append(x_mask[i])

        for i in range(0, len(y_mask)):
            y_mask_ref.append(y_mask[i])

        i = 1
        found = False
        while i < len(y_mask):
            if y_mask[-i] and found:
                y_mask[-i] = False
            elif y_mask[-i]:
                found = True
                y_pos = i
            i += 1

        xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
        self.n[xy_mask] = n_material

        # self._nxt_piece(x_mask, y_mask, n_material, x_mask_ref, y_mask_ref, False)

        return self.n

    def _left_diagonal(self, x_bot_left, y_bot_left, x_top_right, y_top_right,
                     n_material, trap_side = 0, trap_len = 0):

        l_x_max = x_bot_left
        l_x_min = []

        iterations = trap_len/self.x_step
        diff = iterations - int(iterations)

        discrete_step = trap_len / iterations

        trap_range = y_top_right - y_bot_left

        y_max = []
        y_max.append(self.y_min)

        l_x_min = self.x + 1
        new_y = []

        done = False
        i = 1
        while not len(y_max) == len(self.y):
            val = y_max[-1] + (trap_range / iterations)

            while i < len(self.y) and val - self.y[i] > self.x_step and not len(self.y) == len(y_max) + 1:
                # print('comparing:', val, 'against ', self.y[i])
                y_max.append(self.y[i] - 0.0000001)
                i += 1
            if val > y_top_right:
                y_max.append(y_top_right)
            else:
                y_max.append(val)
                new_y.append(i)
                i+=1

        ref_x = self.x
        done = False
        done2 = False
        i = 0

        # print("\nBEFORE:")
        # print('ref_x',ref_x)
        # print('\nl_x_min', l_x_min)


        while i < len(l_x_min):
            if l_x_min[i] - 1 > x_bot_left - trap_len and not done:
                for g in range(0, int(iterations)):
                    l_x_min[i] = (x_bot_left - trap_len) + g * self.x_step
                    i += 1
                    done = True

                if done:
                    l_x_min[i] = x_bot_left
                    if l_x_max % self.x_step != 0:
                        l_x_max += 0.01

                    ref_x[i] = x_bot_left + 0.0001
            i += 1
            #
            # """Questionable"""
            # if (l_x_min[i] - 1 > x_top_right + trap_len and not done2):
            #     for g in range(0, int(iterations)):
            #         l_x_min[i] = (x_top_right + trap_len) - g * self.x_step
            #         i += 1
            #         done2 = True
            #
            #     l_x_min[i] = x_top_right + trap_len
            #     ref_x[i] = x_top_right + trap_len + 0.0001


        # print("\nAFTER:")
        # print('ref_x',ref_x)
        # print('\nl_x_min', l_x_min)

        done = False
        i = 0
        while i < len(y_max) - 1:
            if y_max[i] > self.y[i]:
                while g < int(iterations) and not done:
                    while y_max[i + 1] < self.y[i + 1]:
                        y_max[i + 1] = self.y[i + 1] + 0.01
                        i += 1
                        done = True
            i += 1

        y_top_right = y_max

        x_mask_ref = []
        y_mask_ref = []

        x_mask = np.logical_and(l_x_min<=ref_x, ref_x<=l_x_max)
        y_mask = np.logical_and(y_bot_left<=self.y, self.y<=y_top_right)

        for i in range(0, len(x_mask)):
            x_mask_ref.append(x_mask[i])

        for i in range(0, len(y_mask)):
            y_mask_ref.append(y_mask[i])

        i = 1
        found = False
        while i < len(y_mask):
            if y_mask[-i] and found:
                y_mask[-i] = False
            elif y_mask[-i]:
                found = True
                y_pos = i
            i += 1

        if y_mask[0] and found:
            y_mask[0] = False

        xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
        self.n[xy_mask] = n_material

        self._nxt_piece(x_mask, y_mask, n_material, x_mask_ref, y_mask_ref, True)

        return self.n

    def add_material(self, x_bot_left, y_bot_left, x_top_right, y_top_right,
                     n_material, trap_side = 0, trap_len = 0):

        # print('params')
        print('x_bot_left [x_min]\t', x_bot_left)
        print('y_bot_left \t\t', y_bot_left)
        print('x_top_right \t\t',x_top_right)
        print('y_top_right\t\t', y_top_right)
        # print('\nself params:')
        # print(self.x)
        # print(self.y)

        x_mask = np.logical_and(x_bot_left<=self.x, self.x<=x_top_right)
        y_mask = np.logical_and(y_bot_left<=self.y, self.y<=y_top_right)

        # print('\ny_mask:', y_mask, '\n\n')
        xy_mask = np.kron(y_mask, x_mask).reshape((y_mask.size, x_mask.size))
        self.n[xy_mask] = n_material

        if trap_len:
            self._left_diagonal(x_bot_left, y_bot_left, x_top_right, y_top_right,
                         n_material, trap_side , trap_len)
            self._right_diagonal(x_bot_left, y_bot_left, x_top_right, y_top_right,
                        n_material, trap_side , trap_len)

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

    def add_material(self, x_min, x_max, n, trap_side, trap_len):

        # print('params')
        # print(x_bot_left)
        # print(y_bot_left)
        # print(x_top_right)
        # print(y_top_right)
        # print('self params:')
        # print(self.x)
        # print(self.y)

        """
            Need to make another x_bot_left etc set of params for the diagonal.
            It will need to be a list of values.. Then instead of logical operations, compare in a
            loop such that: x_bot < self.x < x_top. Each must be equal length and aligned.

        """
        # if trap_len and trap_side:
        #     step = trap_len/10
        #     # x_bot[0] = x_bot_left
        #     for i in range (0, trap_len, step):
        #         self.y_min[i] = i/10 * trap_side

        # x_min -= trap_len
        # x_max += trap_len

        # print (self.y_min)
        Structure.add_material(self, x_min, self.y_min, x_max, self.y_max, n, trap_side, trap_len)
        # Structure.left_diagonal(self, x_min, self.y_min, x_max, self.y_max, n, trap_side, trap_len)

        # if trap_len or trap_side:
        #     Structure.add_diagonal_sides(self, x_min, x_max, trap_len, trap_side, n)

        # Structure.add_mat_sides()



        return self.n
