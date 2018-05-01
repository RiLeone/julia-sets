#!/usr/bin/env python3
"""
    Zooming into Julia Sets
    =======================
"""

import juliaSet as js
import numpy as np

if __name__ == "__main__":

    def lin_interp(start, end, nop):
        pts = [start]
        for ii in range(1, nop + 1):
            pts.append(start + ii/nop * (end - start))
        return pts

    print("- Setting up parameters...")
    # Zooming parameters
    start_range_real = 4
    start_range_imag = 4
    start_center = (-0.73214237+0.0001,0.191826425-0.0001)#(0., 0.)
    end_center = (-0.73214237+0.0001,0.191826425-0.0001)
    zooming_factor = 0.95
    zooming_iterations = 120
    start_computation_iterations = 80
    end_computation_iterations = 500

    real_centers = np.asarray(lin_interp(start_center[0], end_center[0], zooming_iterations))
    imag_centers = np.asarray(lin_interp(start_center[1], end_center[1], zooming_iterations))
    comp_iters = np.asarray(lin_interp(start_computation_iterations, end_computation_iterations, zooming_iterations), dtype = 'int')

    real_ranges = [start_range_real]
    imag_ranges = [start_range_imag]
    for ii in range(1, zooming_iterations + 1):
        real_ranges.append(zooming_factor * real_ranges[-1])
        imag_ranges.append(zooming_factor * imag_ranges[-1])

    real_lims = []
    imag_lims = []

    for rerng, imrng, re0, im0 in zip(real_ranges, imag_ranges, real_centers, imag_centers):
        minre = re0 - 0.5 * rerng
        maxre = re0 + 0.5 * rerng

        minim = im0 - 0.5 * imrng
        maxim = im0 + 0.5 * imrng

        real_lims.append((minre, maxre))
        imag_lims.append((minim, maxim))

    C = -0.8 + 0.156 * 1j
    real_res = 1e3
    imag_res = 1e3
    mode = "mandelbrot"

    print("- Starting computation...")
    for idx in range(len(real_lims)):
        print("\tComputing iteration {:d} of {:d}...".format(idx + 1, len(real_lims)))
        framename = "frames/{:05d}.png".format(idx + start_computation_iterations)
        julia = js.juliaSet(C, real_res, imag_res, real_lims[idx], imag_lims[idx], mode, False)
        if idx == 0:
            julia.compute_iterations(comp_iters[idx], True)
        else:
            julia.compute_iterations(comp_iters[idx], False)
            julia.draw_result(savefig = True, savefilename = framename, showdrawing = False)

        del(julia)


    print("* ============= *")
    print("* End of Script *")
    print("* ============= *")
