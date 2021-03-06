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
import PIL.Image

class juliaSet:
    def __init__(self,
                 C = -1.,
                 real_res = 1e2,
                 imag_res = 1e2,
                 real_lim = (-1.1, 1.1),
                 imag_lim = (-1.1, 1.1),
                 mode = "quadratic",
                 verbose = True):

        self._real_res = real_res
        self._imag_res = imag_res
        self._real_lim = real_lim
        self._imag_lim = imag_lim
        self._mode = mode
        self._verbose = verbose
        self._noi = 0
        self._real_range = self._real_lim[1] - self._real_lim[0]
        self._imag_range = self._imag_lim[1] - self._imag_lim[0]
        self._range_ratio = self._imag_range / self._real_range

        self._plane_mesh = np.meshgrid(np.linspace(self._real_lim[0], self._real_lim[1], self._real_res), np.linspace(self._imag_lim[0], self._imag_lim[1], self._imag_res))

        self._z = self._plane_mesh[0] + 1j * self._plane_mesh[1]

        if mode == "mandelbrot":
            self._C = self._z
        else:
            self._C = C * np.ones_like(self._z)

        self._julia = np.float64(np.zeros_like(np.real(self._z)))

    def compute_iterations(self, noi = 1e2, drawandsaveframes = False, filename = None):

        self._noi = noi

        for ii in range(int(noi)):
            if self._verbose:
                print("Computing iteration {:d} of {:d}".format(ii + 1, noi))

            if drawandsaveframes:
                if filename == None:
                    framefilename = "frames/{:05d}.png".format(ii)
                else:
                    framefilename = filename

                self.draw_result(visualization = "imshow", savefig = True, savefilename = framefilename, showdrawing = False)

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

        # if drawandsaveframes:
        #     if filename == None:
        #         framefilename = "frames/{:05d}.png".format(noi)
        #     else:
        #         framefilename = filename
        #
        #     self.draw_result(visualization = "imshow", savefig = True, savefilename = framefilename, showdrawing = False)


    def draw_result(self,
                    colormap = cm.hot,
                    visualization = "imshow",
                    savefig = False,
                    savefilename = None,
                    showdrawing = True):

        if not visualization == "PIL":
            fig = pltlib.figure(figsize = (int(self._real_res / 100), self._range_ratio * int(self._real_res / 100)))
            ax = fig.add_subplot(111, aspect = 'equal')
            ax.set_position([0., 0., 1., 1.])

        if savefilename == None:
            filename = visualization + "_"
            filename += "julia_set_mode_" + self._mode
            if not mode == "mandelbrot":
                filename += "_C_{:.2f}j{:.2f}".format(np.real(self._C[0][0]), np.imag(self._C[0][0]))
            filename += "_re{:.2e}_{:.2e}".format(self._real_lim[0], self._real_lim[1])
            filename += "_im{:.2e}_{:.2e}".format(self._imag_lim[0], self._imag_lim[1])
            filename += "_iters{:d}.png".format(self._noi)
        else:
            filename = savefilename

        if visualization == "imshow":
            pltlib.imshow(self._julia,
                          cmap = colormap,
                          interpolation = 'bicubic',
                          aspect = 'equal',
                          origin = 'lower',
                          vmin = 0.,
                          vmax = 1.,
                          extent = (0., 1., 0., 1.))
            pltlib.axis('off')

        elif visualization == "PIL":
            a = self._julia
            a_cyclic = (6.28 * a / 20.).reshape(list(a.shape) + [1])
            # img = np.concatenate([ 10 + 20 * np.cos(a_cyclic),
            #                        30 + 50 * np.sin(a_cyclic),
            #                       155 - 80 * np.cos(a_cyclic)], 2)
            img = np.concatenate([50 - 50 * np.sin(a_cyclic),
                                  10 - 10 * np.sin(a_cyclic),
                                  155 - 100 * np.sin(a_cyclic)], 2)
            # img[a == a.max()] = 255
            a = img
            a = np.uint8(np.clip(a, 0, 255))
            img = PIL.Image.fromarray(a)
            img.save("plts/" + filename)
            img.close()# img.show()

        elif visualization == "surf":
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
            if not visualization == "PIL":
                pltlib.savefig("plts/" + filename)
            else:
                pass # img.save('plts/' + filename)

        if showdrawing:
            if not visualization == "PIL":
                pltlib.show()
            else:
                pass # img.show()
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

    phi = (1. + np.sqrt(5)) / 2.
    C = -0.8 + 0.156 * 1j
    real_res = 2e3
    imag_res = 2e3
    # Explorer parameters
    x0 = -0.8
    y0 = 0.1
    span0 = 0.1

    center_factors = [[0.627083, 0.829167], [0.17778, 0.546528], [0.674306, 0.688889]]
    new_spans = [0.01, 0.005, 0.002]
    for cf, ns in zip(center_factors, new_spans):
        x0 += cf[0] * span0
        y0 += cf[1] * span0
        span0 = ns

    print("Center: (x0, y0) = ({},{})\nSpan: span0 = {}".format(x0, y0, span0))

    real_lim = (x0 - 0.5 * span0, x0 + 0.5 * span0) # (-0.1, 0.1)# 1.7  #(-0.8, -0.7)
    imag_lim = (y0 - 0.5 * span0, y0 + 0.5 * span0)  # (-0.1, 0.1)# 1.2  #( 0.1, 0.2)

    real_lim = (-1.7, 1.7)  #(-0.8, -0.7)
    imag_lim = (-1.2, 1.2)  #( 0.1, 0.2)

    mode = "quadratic"

    print("- Initializing object...")
    test = juliaSet(C = C,
                    real_res = real_res,
                    imag_res = imag_res,
                    real_lim = real_lim,
                    imag_lim = imag_lim,
                    mode = mode)

    print("- Start computing iterations...")
    test.compute_iterations(200, drawandsaveframes = False)

    print("- Showing final frame...")
    test.draw_result(visualization = "imshow", savefig = True)

    print("* ============= *")
    print("* End of Script *")
    print("* ============= *")
