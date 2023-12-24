# Bispectrum

This repository houses codes used by the K-State cosmology group to compute
bispectrum of galaxies from various simulations and data.

The release 1.0 computes bispectrum multipoles from uniform periodic cubic boxes.
The key functions are tested.

The main function for end users is in "scripts/compute_bispectrum.jl"
compute_pk_bk(x, y, z, dk, N, L, outfile)
x, y, z are the vectors with coordinates. They must be shifted so that they all
start from 0!!!
dk is the bin width for the bispectrum multipoles in k1, k2, k3
N is the number of bins
L is the size of the box (in Mpc/h)
The bispectrum will be written in a file "bk_" + outfile
The file will have 6 columns for k1, k2, k3, B0, B2, B4
The powerspectrum will be written in a file "pk_" + outfile
The file will ahve 4 columns for k, P0, P2, P4

The code is multithreaded and can be run with e.g. "julia --threads 128 xxx.jl"

The performance scales very steeply with dk*N.
For dk = 0.01 and N = 30 it takes about 10 minutes on a single Perlmutter node
with 128 threads (hopefully this will improve significantly in the next
release!)

Things that will be done by the next release:

* Computing Bispectrum from non-uniform surveys (using random catalogs)
* Accounting for the wide-angle effects in the multipoles
* Full coverage of all functions
* A potentially much faster algorithm based on recent ideas
