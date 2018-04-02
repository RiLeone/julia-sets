# Julia Sets - Nerd Art with Python

Computation and visualization of Julia-sets with Python. There are probably a million ways
of doing this more efficiently, but this version works.

Check [Julia Set on Wiki](https://en.wikipedia.org/wiki/Julia_set) for an intro on these beauties and for understanding why the Mandelbrot set can be represented with this library as well.

## Use

The Julia set is implemented as a custom class in the file `juliaSet.py`, this implements a class `juliaSet`, which has methods for computation and visualization of the set.

The file `zooming.py` uses instances of the `juliaSet` object to create frames for animations zooming into a specific region of a specific set.

## Remark on dummy files

In order to be able to ignore `.png`s and `.gif`s while having the needed folder structure,
`dummy.md` files are placed in the required subfolders.
