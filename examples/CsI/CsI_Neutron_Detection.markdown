# Neutron Detection Approaches in CsI

In this example, we aim to explore approaches suitable for neutron detection in CsI. In CsI, there are only two atoms per unit cell, and the high-symmetry-direction branches are purely longitudinal or transverse.

Below, we show the process of calculating the dynamical structure factor (SQE) for a unit cell with two atoms.

### Input Parameters Setup

For example, we want to explore longitudinal acoustic phonon along [H 0 0],thus we condisering as below:

### 1. Basic Section
Based on the POSCAR entry `"Cs I 1 1"`, we define:

- **ntypes = 2**: Two atom types (Cs and I).
- **natoms = 2**: Two atoms in the unit cell.
- **nsize = 2 2 2**: Supercell size for FORCE_CONSTANTS calculation is 2x2x2.

Since CsI has a cubic lattice with lattice vectors:  
```
4.5199999809 0.0000000000 0.0000000000
0.0000000000 4.5199999809 0.0000000000
0.0000000000 0.0000000000 4.5199999809
```
the `clatvec` parameter can be omitted in this case.

### 2. Input SQE Section
From the POSCAR entry `"Cs I 1 1"`, we set:  
- **elements = "Cs" "I"**: Specifies the elements.
- **nat = 1 1**: One atom of each type.

For neutron scattering calculations:  
- **lneutron = .TRUE.**: Enable neutron scattering.
- **lxray = .FALSE.**: Disable X-ray scattering.

To ensure converged root-mean-square displacement (RMSD):  
- **qmesh = 20 20 20**: Q-point mesh for integration.
- **write_rmsd = .TRUE.**: Write RMSD in the first calculation.
- **read_rmsd = .TRUE.**: For subsequent calculations, read RMSD (commented out in the first run).

To calculate along the [H 0 0] direction:  
- **path(:,1) = 1 0 0**: Specifies the [H 0 0] direction. The other directions (`path(:,2)` and `path(:,3)`) are automatically determined.  
  Note: For CsI’s cubic lattice, [1 0 0] corresponds to [H 0 0].

For energy settings with a 12 meV cutoff:  
- **ne = 240**: Number of energy points.
- **deltaE = 0.05**: Energy step size (meV).

For SQE calculation along the H direction:  
- **nqh = 400**: Number of Q-points along H.
- **deltaH = 0.01**: Step size along H.

For integration in the K and L directions:  
- **nqk = 5**, **nql = 5**: Number of Q-points in K and L directions.
- **deltaK = 0.005**, **deltaL = 0.005**: Step sizes in K and L directions.  
  Total integration thickness:  
  - K direction: `nqk * deltaK = 5 * 0.005 = 0.025`  
  - L direction: `nql * deltaL = 5 * 0.005 = 0.025`

Starting point for the Q-vector:  
- **q0 = 0 0 0**: Origin of the reciprocal lattice.

No Born correction is considered:  
- **nonanalytic = .FALSE.**

Instrument resolution for CNCS20:  
- **lresfunc = .TRUE.**: Enable resolution function.
- **functype = "CNCS20"**: Specify CNCS20 instrument.
- **xm = 1**: Instrument parameter.

### 3. Lattice Parameters
Based on the POSCAR, the lattice parameters are defined as:  
```
LATTICE_PARAMETERS
1.d0
4.5199999809 0.0000000000 0.0000000000
0.0000000000 4.5199999809 0.0000000000
0.0000000000 0.0000000000 4.5199999809
```

### 4. Atomic Positions
The atomic positions are:  
```
ATOMIC_POSITIONS
Direct
0.000000000 0.000000000 0.000000000
0.500000000 0.500000000 0.500000000
```

### Complete Input File
Below is the complete input file for the SQE calculation:

```
! SQE input file for CsI

&basic
ntypes = 2 
natoms = 2
nsize = 2 2 2
/

&inputsqe
elements = "Cs" "I"
nat = 1 1 
temp = 0.d0
lneutron = .TRUE.
lxray = .FALSE.
qmesh = 20 20 20
write_rmsd = .TRUE.
!read_rmsd = .TRUE.
path(:,1) = 1 0 0
ne = 240
deltaE = 0.05
nqh = 400
nqk = 5
nql = 5
deltaH = 0.01
deltaK = 0.005
deltaL = 0.005
q0 = 0 0 0
nonanalytic = .FALSE.
lresfunc = .TRUE.
functype = "CNCS20"
xm = 1
/

LATTICE_PARAMETERS
1.d0
4.5199999809 0.0000000000 0.0000000000
0.0000000000 4.5199999809 0.0000000000
0.0000000000 0.0000000000 4.5199999809

ATOMIC_POSITIONS
Direct
0.000000000 0.000000000 0.000000000
0.500000000 0.500000000 0.500000000
```

## Execution and Visualization
1. Execute the command `dysurf CsI.txt` to compute the SQE data.
2. Use the Python script `h00plot.py`(at "./h00/") to plot the SQE results along the [H 0 0] direction.

This setup enables the calculation and visualization of the dynamical structure factor for neutron detection in CsI, focusing on the [H 0 0] direction with appropriate instrument resolution to explore the longitudinal acoustic phonon modes.

The result is shown below: 
![SQE along [H 0 0]](./h00/SQE_0K.png)

## Get SQE of other direction:
The results for the SQE calculation are shown in the figure below:
To obtain SQE results along different paths, we primarily adjust path(:,1) and q0. For example:

Setting path(:,1) = 1 0 0 and q0 = 0 1 0 allows us to compute the transverse acoustic phonon along the [H 0 0] direction. The input file is located at ./h10/Cs.txt. The result is shown in the figure below:
![SQE along [H 1 0]](./h10/SQE_0K.png)
Setting path(:,1) = 1 0 0 and q0 = 0 2 0 allows us to compute the transverse acoustic phonon along the [H 0 0] direction. The input file is located at ./h20/Cs.txt. The result is shown in the figure below:
![SQE along [H 2 0]](./h20/SQE_0K.png)
