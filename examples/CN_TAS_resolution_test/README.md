# CN_TAS_resolution_test

This directory collects three TAS resolution examples for Dysurf, which relies on Cooper–Nathans method. They are arranged from the most local diagnostic case to the S(**Q**,E) comparison workflow.

## Subdirectories

- `2+H00_psf2d_local_example`
  A minimal local PSF example at one representative point. It is mainly used to check that the local TAS 4D resolution can be approximated by the fitted 2D PSF form.

- `2+H00_psf2d_hcol_only_compare`
  A local PSF comparison where only the horizontal collimations are changed. It is used to show how `hcol` changes the fitted PSF shape while keeping the other TAS parameters fixed.

- `2+H00_CNCS12_psf2d_hcol_compare`
  A S(**Q**,E) comparison example. It compares the legacy empirical `CNCS12` resolution path with two formal TAS `psf2d` calculations using different horizontal collimations.

## Suggested Reading Order

If you want to understand the workflow step by step, the recommended order is:

1. `2+H00_psf2d_local_example`
2. `2+H00_psf2d_hcol_only_compare`
3. `2+H00_CNCS12_psf2d_hcol_compare`

The first two folders focus on local resolution diagnostics. The third folder shows how those ideas appear in a full `S(Q,E)` calculation.
