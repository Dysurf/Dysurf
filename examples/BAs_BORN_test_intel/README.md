# BAs_BORN_test

This folder is a minimal standalone smoke test for the `BORN` / `nonanalytic` branch in the Intel build of Dysurf.

It reuses the same BAs Born-enabled input physics as the main [BAs](/home/davy/software/Dysurf/examples/BAs) example, but keeps only one input file so it can be used as a quick regression test for `src/src_intel`.

## Files

- [BAs_BORN_test.txt]
- [FORCE_CONSTANTS]
- [run_example.sh]

## What It Tests

- `nonanalytic = .TRUE.`
- dielectric tensor from the `BORN` block
- Born effective charges from the `BORN` block
- BAs phonon and SQE generation with the current [src/src_intel](/home/davy/software/Dysurf/src/src_intel) build

## How To Run

```bash
cd Dysurf/examples/BAs_BORN_test
bash run_example.sh
```

This writes the standard Dysurf outputs in the current folder, including:

- `SQE_300K.dat`
- `omega.dat`
- `qpoints.dat`
- `qpoints_full.dat`
- `rmsd.dat`

## Notes

- [FORCE_CONSTANTS](/home/davy/software/Dysurf/examples/BAs_BORN_test/FORCE_CONSTANTS) is a symbolic link to the shared BAs force-constant file in [examples/BAs](/home/davy/software/Dysurf/examples/BAs).
- This folder is intended as a compact Intel-build functional test, not as the full Born vs no-Born comparison workflow.
