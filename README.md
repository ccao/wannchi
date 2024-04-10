# WannChi: Susceptibility from Wannier Hamiltonians

This program was developed by Chao Cao, Siqi Wu, Chenchao Xu, and Guo-Xiang Zhi @ Zhejiang University. 

Bugs and comments can be submitted via GitHub or email to:

ccao@zju.edu.cn

## Installation

Some implementation of BLAS/LAPACK is required, MKL is tested. MPI is recommended.

1. Modify "make.sys" according to your specification.

2. Get into modules directory and make.

3. Get into src directory and make.

Repeat 1-3 if you encountered compilation problems. 

## Run the program

Prepare some input files and use:

```
    mpirun -np ${NUM_PROC} wannchi.x > wannchi.log
```

## Input Files:

* wannchi.inp

    Main input file.

* _seed_\_hr.dat

    Wannier Hamiltonian from Wannier90.

* _seed_.pos

    POSCAR-like input for atomic positions.

* IBZKPT

    VASP KPOINTS-like input for K-mesh definition.

* QPOINTS

    Specifies q-points to be calculated.

* RPA.inp

    Optional, if full matrix calculation is required.

### Description of _seed_.pos file

POSCAR-like input for orbital-site information.

Format description:
```
line 1            : Some comment
line 2            : Universal Scale
line 4            : A1
line 5            : A2
line 6            : A3
line 7            :  nsite  soc
line 8 to nsite+7 : Zat  x1  x2  x3   nbasis
```
nsite: Number of sites
soc  : 0 (no SOC) or 1 (with SOC)
Zat  : Atomic number
x1, x2, x3 : fractional coordinates
nbasis : Number of Wannier orbitals on this site

Example 1: Perfect Kagome w/o SOC
```
Kagome
 1.0
     4.6669313530311989   -2.6944540729627522   -0.0000000000000000
     4.6669313530311989    2.6944540729627526   -0.0000000000000000
     0.0000000000000000    0.0000000000000000    9.8872122322434954
  3  0
1  0.5  0.0  0.0 1
1  0.0  0.5  0.0 1
1  0.5  0.5  0.0 1
```

### Description of IBZKPT file

KPOINTS-like input for K-mesh.

Format description:
```
line 1 : Some comment
line 2 : switch
line 3 : Some other comment
  nk1  nk2  nk3
```
switch :   0 for automatic generation. Currently only 0 is supported.
nk1, nk2, nk3 : K-mesh definition. Currently only Gamma-centered K-mesh is supported.

Example : A 96 x 96 x 1 Gamma-centered K-mesh
```
Automatically generated mesh
      0
Reciprocal lattice
  96  96  1
```

### Description of QPOINTS file

Specifies q-points to be calculated. The format is controlled by the first line, which is a mode switch. Possible values are:

* 0 : single-point calculation. The 2nd line specifies the q-point. All other lines are ignored.
* 1 : line-mode calculation. This is useful to generate band structure like data.
* 2 : plane-mode calculation. Calculate on a 2D mesh.
* 3 : bulk-mode calculation. 

Example 1 : Single point calculation on (\pi, \pi, 0)

```
0
0.5  0.5  0.0
```

Example 2 : Line mode calculation for square lattice along G-M-X-G. Each segment is interpolated by 24 points.
```
1
3  24
0.0  0.0  0.0   0.5  0.5  0.0
0.5  0.5  0.0   0.5  0.0  0.0
0.5  0.0  0.0   0.0  0.0  0.0
```

Example 3 : Plane-mode calculation for a 48x48 mesh spanned by vectors (1.0, 0.0, 0.0)  and (0.0, 1.0, 0.0), with one end at (0.0, 0.0, 0.5)
```
2
48 48
0.0  0.0  0.5
1.0  0.0  0.0
0.0  1.0  0.0
```

Example 4 : Body-mode calculation for a 24x24x24 mesh.
```
3
24 24 24
```

### Description of RPA.inp file

### Description of wannchi.inp
