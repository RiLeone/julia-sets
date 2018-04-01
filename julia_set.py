#!/usr/bin/env python3
"""
    Julia Sets
    ==========
"""
import numpy as np
import matplotlib.pyplot as pltlib
import matplotlib.colors as clrs
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class juliaSet:
    def __init__(self,
                 C = -1.,
                 real_res = 1e2,
                 imag_res = 1e2,
                 real_lim = (-1.1, 1.1),
                 imag_lim = (-1.1, 1.1),
                 mode = "quadratic"):

        self._real_res = real_res
        self._imag_res = imag_res
        self._real_lim = real_lim
        self._imag_lim = imag_lim
        self._mode = mode
        self._noi = 0

        self._plane_mesh = np.meshgrid(np.linspace(self._real_lim[0], self._real_lim[1], self._real_res), np.linspace(self._imag_lim[0], self._imag_lim[1], self._imag_res))

        self._z = self._plane_mesh[0] + 1j * self._plane_mesh[1]

        if mode == "mandelbrot":
            self._C = self._z
        else:
            self._C = C * np.ones_like(self._z)

        self._julia = np.float64(np.zeros_like(np.real(self._z)))

    def compute_iterations(self, noi = 1e2, drawandsaveframes = False):

        self._noi = noi

        for ii in range(int(noi)):
            print("Computing iteration {:d} of {:d}".format(ii + 1, noi))

            if drawandsaveframes:
                framefilename = "frames/{:05d}.png".format(ii)
                self.draw_result(useimshow = True, savefig = True, savefilename = framefilename, showdrawing = False)

            if self._mode == "exponential":
                self._z = np.exp(self._z) + self._C
            elif self._mode == "quadratic":
                self._z = self._z ** 2 + self._C
            elif self._mode == "quadlog":
                self._z = (self._z ** 2 + self._z) / np.log(self._z) + self._C
            else:
                self._z = self._z ** 2 + self._C

            self._julia += np.float64((np.real((self._z * np.conjugate(self._z))) > 4) / (ii + 1))

        self._z = np.nan_to_num(self._z)
        self._julia = np.nan_to_num(self._julia)

        if drawandsaveframes:
            framefilename = "frames/{:05d}.png".format(noi)
            self.draw_result(useimshow = True, savefig = True, savefilename = framefilename, showdrawing = False)


    def draw_result(self, colormap = cm.hot, useimshow = True, savefig = False, savefilename = None, showdrawing = True):

        fig = pltlib.figure(figsize = (6,6))

        if useimshow:
            pltlib.imshow(self._julia,
                          cmap = colormap,
                          interpolation = 'bicubic',
                          aspect = 'equal',
                          origin = 'lower',
                          vmin = 0.,
                          vmax = 1.,
                          extent = (0., 1., 0., 1.))
            pltlib.axis('off')
        else:
            ax = fig.add_subplot(111,
                                 projection='3d',
                                 aspect = 'equal')

            ax.plot_surface(self._plane_mesh[0],
                            self._plane_mesh[1],
                            self._julia,
                            ccount = 512,
                            cmap =  colormap,
                            shade = True,
                            linewidth = 0,
                            antialiased = False)

            # ax.set_xlabel("real")
            ax.view_init(90, -90)
            pltlib.subplots_adjust(top = 1,
                                   bottom = 0.,
                                   left = 0.0,
                                   right = 1.,
                                   hspace = 0.,
                                   wspace = 0.)
            ax.axis('off')

        if savefig:
            if savefilename == None:
                filename = "julia_set_mode_" + self._mode
                if not mode == "mandelbrot":
                    filename += "_C_{:.2f}j{:.2f}".format(np.real(self._C[0][0]), np.imag(self._C[0][0]))
                filename += "_re{:.2e}_{:.2e}".format(self._real_lim[0], self._real_lim[1])
                filename += "_im{:.2e}_{:.2e}".format(self._imag_lim[0], self._imag_lim[1])
                filename += "_iters{:d}.png".format(self._noi)
            else:
                filename = savefilename

            pltlib.savefig("plts/" + filename)

        if showdrawing:
            pltlib.show()
        else:
            pltlib.close(fig)


if __name__ == "__main__":

    print("Julia Sets Visualization")
    print("========================\n")
    # Interesting values of C when not going for the Madelbrot option:
    #   -0.8 + 0.156 * 1j
    #   0.65 * 1j
    #   1 - phi = 1 - (1 + np.sqrt(5)) / 2
    #   phi - 2 + 1j * (phi - 1)
    #   0.285 + 0.01 * 1j
    # For exp version
    #   -0.65
    # For quadlog version
    #   0.268 + 0.060 * 1j

    C = 0.268 + 0.060 * 1j
    real_res = 1e3
    imag_res = 1e3
    real_lim = (-1.7, 1.7)
    imag_lim = (-1.2, 1.2)
    mode = "quadlog"

    print("- Initializing object...")
    test = juliaSet(C = C,
                    real_res = real_res,
                    imag_res = imag_res,
                    real_lim = real_lim,
                    imag_lim = imag_lim,
                    mode = mode)

    print("- Start computing iterations...")
    test.compute_iterations(100, drawandsaveframes = True)

    print("- Showing final frame...")
    test.draw_result(useimshow = True, savefig = True)

    print("* ============= *")
    print("* End of Script *")
    print("* ============= *")
