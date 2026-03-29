# background
The measured TAS intensity is affected not only by the intrinsic 
𝑆(𝑄, 𝐸), but also by the instrumental resolution. In this code, two related descriptions are provided.

The Cooper–Nathans (CN) method is the standard instrument-based TAS resolution formalism. It is useful for understanding how the experimental geometry, collimations, mosaics, and fixed-energy settings determine the resolution.

The local 2D PSF method is a reduced representation of the local TAS resolution in the slice plane. It is useful for efficient convolution, visualization, and comparison between intrinsic and broadened 
𝑆(𝑄, 𝐸) slices.

In short, CN is the more physical instrument-level description, while the PSF method is the more practical slice-level representation.

# 2+H00_psf2d_local_example

This folder contains a single local `psf2d` validation example for the PdSe2 `2+H00` scan. Its purpose is to show, at one representative `(Q, E)` point, how the raw local TAS Cooper-Nathans resolution patch compares with the fitted sheared 2D PSF used by Dysurf.

In other words, this is the smallest example for checking whether the `psf2d` parameterization can reproduce the local TAS response shape.

## What This Folder Does

- chooses one representative point automatically from the base shear field
- exports the raw local 2D PSF around that point
- exports the fitted PSF built from `sigma_q`, `sigma_e_left`, `sigma_e_right`, `shear`
- plots raw map, fitted map, `dQ` cut, and `dE` cut

## TAS Parameters

The TAS-related parameters are documented in [instrument_base.txt](/home/davy/software/Dysurf/examples/CN_TAS_resolution_test/2+H00_psf2d_local_example/instrument_base.txt) with comments. The main ones are:

- `tas_mode = 'CN'`: use Cooper-Nathans TAS resolution
- `tas_obs_mode = 'psf2d'`: export the fitted local 2D PSF representation
- `tas_fix_mode = 'Ef'`: fixed final energy mode
- `tas_e_fixed = 14.7d0`: final neutron energy in meV
- `tas_coll_h_*`: horizontal collimations in arcmin
- `tas_mosaic_*_h`: horizontal mosaics in arcmin
- `tas_dm`, `tas_da`: monochromator / analyzer d spacing
- `u`, `v`: sample orientation vectors

`u` and `v` are now real runtime inputs. `run_psf2d_example.py` reads them from [instrument_base.txt](/home/davy/software/Dysurf/examples/CN_TAS_resolution_test/2+H00_psf2d_local_example/instrument_base.txt) and passes them to [validate_psf2d_parametrization.f90](/home/davy/software/Dysurf/src/src_gfortran/validate_psf2d_parametrization.f90).
Their meaning is now exactly the same as the `u/v` parameters in the SQE `psf2d` example: they define the sample orientation used by the local TAS 4D-resolution to PSF parameter pipeline, not the scan path itself.

For the current example, the orientation is:

- `u = [1, 0, 0]`
- `v = [0, 1, 0]`

## How To Run

```bash
cd Dysurf/examples/CN_TAS_resolution_test/2+H00_psf2d_local_example
bash run_example.sh
```

## Figure

The upper panels compare the original TAS resolution function obtained from the Cooper–Nathans (CN) method with the corresponding fitted sheared 2D PSF. Both show a very similar tilted elliptical shape in the local (𝑑𝑄, 𝑑𝐸) plane, indicating that the fitted PSF captures the main features of the original TAS resolution, including the widths in momentum and energy as well as the 𝑄–𝐸 correlation (shear/tilt).

The lower panels show one-dimensional cuts along the momentum and energy directions. The fitted curves agree closely with the CN results, demonstrating that the local 2D PSF provides a good reduced representation of the full CN resolution for practical slice-based convolution and visualization.

In short, the CN method gives the more complete instrument-based resolution description, while the fitted local 2D PSF offers a simpler and efficient approximation that is convenient for 2D 
𝑆(𝑄, 𝐸) calculations.

- [psf2d_example.png](./psf2d_example.png)