# Tutorial: SQE Calculation for GaN

This tutorial outlines the process of setting up the input parameters for calculating the dynamical structure factor (SQE) for GaN along specific high-symmetry directions, focusing on slice line profiles along [H 0 0]. The input file provided is tailored for GaN with a wurtzite structure (space group 186), and we explain each parameter section in detail.

## Parameter Settings

### 1. Basic Section

The `&basic` section defines the fundamental properties of the system based on the structural parameters in the POSCAR file:

- **ntypes = 2**: Specifies two atom types, Ga (Gallium) and N (Nitrogen).
- **natoms = 4**: Indicates four atoms in the unit cell, consistent with the wurtzite structure of GaN.
- **nsize = 4 4 3**: Defines a 4x4x3 supercell for the FORCE_CONSTANTS calculation, which ensures sufficient sampling of phonon interactions.

### 2. Input SQE Section

The `&inputsqe` section configures the SQE calculation parameters:

- **elements = "Ga" "N"**: Specifies the elements in the unit cell (Gallium and Nitrogen).
- **nat = 2 2**: Indicates two atoms of each type (Ga and N) in the unit cell.

To perform the calculation in an orthogonal lattice, we define the conventional lattice vectors (`clatvec`) as:

```
clatvec(:,1) = 3.1900000572   0.0000000000   0.0000000000
clatvec(:,2) = 0.0000000000   5.5252421752   0.0000000000
clatvec(:,3) = 0.0000000000   0.0000000000   5.1900000572
```

These vectors are derived from the primitive cell lattice vectors in the POSCAR:

```
3.1900000572   0.0000000000   0.0000000000
-1.5950000286   2.7626210876   0.0000000000
0.0000000000   0.0000000000   5.1900000572
```

- **temp = 300.d0**: Sets the temperature to 300 K for the calculation, reflecting typical experimental conditions.
- **lneutron = .TRUE.**: Enables neutron scattering calculations.
- **lxray = .FALSE.**: Disables X-ray scattering calculations.
- **qmesh = 20 20 20**: Specifies a 20x20x20 Q-point mesh for integration to ensure converged root-mean-square displacement (RMSD).
- **write_rmsd = .TRUE.**: Writes the RMSD data in the first calculation. For subsequent runs, you can comment this out and use `read_rmsd = .TRUE.` to read the previously computed RMSD.

The high-symmetry path is defined with:

- **path(:,1) = 1 1 0**: Specifies the primary direction for the SQE calculation, which corresponds to a high-symmetry direction in the crystal structure. In this case, the direction is given as [H 0 0] in the conventional cell. For guidance on converting directions defined under LATTICE_PARAMETERS to the clatvec format, please refer to the tutorial located in "../3D_SQE/".
- **path(:,2) = 0 2 0**: Defines the second direction for the path (optional, used for additional sampling or validation).
- **path(:,3) = 0 0 1**: Defines the third direction, typically along the c-axis in wurtzite GaN.

path(:,1), path(:,2) and path(:,3) are mutually perpendicular.

Energy settings:

- **ne = 100**: Sets 100 energy points for the calculation.
- **deltaE = 0.5**: Specifies an energy step size of 0.5 meV, covering an energy range of 50 meV (100 * 0.5).

Q-point sampling along the H direction:

- **nqh = 100**: Number of Q-points along the H direction.
- **deltaH = 0.005**: Step size along H, resulting in a total propagation distance of `nqh * deltaH = 100 * 0.005 = 0.5`.

Integration in the K and L directions:

- **nqk = 10**, **nql = 10**: Number of Q-points in the K and L directions for integration.
- **deltaK = 0.01**, **deltaL = 0.01**: Step sizes in K and L directions, yielding integration thicknesses of `nqk * deltaK = 10 * 0.01 = 0.1` and `nql * deltaL = 10 * 0.01 = 0.1`.

Starting point for the Q-vector:

- **q0 = 0 0 4**: Sets the origin of the reciprocal lattice path.

No LO-TO splitting is considered:

- **nonanalytic = .FALSE.**

Instrument resolution settings for CNCS12:

- **lresfunc = .TRUE.**: Enables the instrument resolution function.
- **functype = "CNCS12"**: Specifies the CNCS12 instrument.
- **xm = 2**: Sets the instrument parameter.

Phase factor inclusion:

- **lphase = .TRUE.**: Enables phase factor calculations, which may be relevant for capturing interference effects in neutron scattering. Here, there is little effect on the results of GaN

### 3. Lattice Parameters

The lattice parameters are based on the POSCAR and define the primitive cell:

```
LATTICE_PARAMETERS
1.d0
3.1900000572   0.0000000000   0.0000000000
-1.5950000286   2.7626210876   0.0000000000
0.0000000000   0.0000000000   5.1900000572
```

### 4. Atomic Positions

The atomic positions in direct coordinates, as specified in the POSCAR, are:

```
ATOMIC_POSITIONS
Direct
0.333333357   0.666666714   0.000000000
0.666666620   0.333333314   0.500000000
0.333333357   0.666666714   0.376500004
0.666666620   0.333333314   0.876500004
```

### Complete Input File

Below is the complete input file for the SQE calculation:

```
! Example of SQE input file for GaN

&basic
ntypes = 2 
natoms = 4
nsize = 4 4 3
/

&inputsqe
elements = "Ga" "N"
nat = 2 2
clatvec(:,1) = 3.1900000572   0.0000000000   0.0000000000
clatvec(:,2) = 0.0000000000   5.5252421752   0.0000000000
clatvec(:,3) = 0.0000000000   0.0000000000   5.1900000572 
temp = 300.d0
lneutron = .TRUE.
lxray = .FALSE.
qmesh = 20 20 20
write_rmsd = .TRUE.
!read_rmsd = .TRUE.
path(:,1) = 1 1 0
path(:,2) = 0 2 0
path(:,3) = 0 0 1
ne = 100
deltaE = 0.5
nqh = 100
nqk = 10
nql = 10
deltaH = 0.005
deltaK = 0.01
deltaL = 0.01
q0 = 0 0 4
nonanalytic = .FALSE.
lresfunc = .TRUE.
functype = "CNCS12"
xm = 2
lphase = .TRUE.
/

LATTICE_PARAMETERS
1.d0
3.1900000572   0.0000000000   0.0000000000
-1.5950000286   2.7626210876   0.0000000000
0.0000000000   0.0000000000   5.1900000572

ATOMIC_POSITIONS
Direct
0.333333357   0.666666714   0.000000000
0.666666620   0.333333314   0.500000000
0.333333357   0.666666714   0.376500004
0.666666620   0.333333314   0.876500004
```

### Execution

1. Place the input file (e.g., `GaN.txt`) in the working directory.
2. Run the command:
   ```bash
   /your_path_to_dysurf/dysurf GaN.txt
   ```
   This generates the SQE data based on the specified parameters.

### Extract the line profile.
We provide a script "./h04/plot_lines.py" to extract line profile at different q points. We just need to adjust the integrated bin range at around line 18 in "./h04/plot_lines.py". And then execute 
```bash
cd h04
python plot_lines.py
```
we will get result as shown below:
![SQE at 300 K along [H 0 0]](./h04/TA2_in_GaN_along_H04.png)

This setup enables the calculation of the dynamical structure factor and get line profile for GaN, capturing phonon modes along the specified [H 0 0] direction with neutron scattering and appropriate instrument resolution.