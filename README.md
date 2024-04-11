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

Example 1 : Single point calculation on ($\pi$, $\pi$, 0)

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

Required if trace_only tag is set to .false.. It specifies the "FF" block as well as the interactions.

Format description:
```
line 1 : nFFblk : Number of FF blocks
line 2 : Dimension of the 1st FF block
line 3 : Global index of the 1st FF block orbitals
line 4 : Dimension of the 2nd FF block
line 5 : Global index of the 2nd FF block orbitals
...
line 2*nFFblk : Dimension of the last FF block
line 2*nFFblk+1 : Global index of the last FF block orbitals
line 2*nFFblk+2 : nUint : Number of interactions specified, can be 0
line 2*nFFblk+3 to line 2*nFFblk+nUint+2 : i1 i2  j1 j2  Uint
```

RPA feature is currently not released to public. The meaning of FF-trick will be explained in future paper.

Example 1: To obtain full matrix for Perfect Kagome.
```
1
3
1  2  3
0
```

Example 2: Full matrix calculation for iron-pnictides with first 10 orbitals being Fe-3d
```
2
 5
1 2 3 4 5
 5
6 7 8 9 10
0
```

### Description of wannchi.inp

Main input file. Fortran namelist input format. Requires 2 blocks of namelist, namely "system" and "control".

Typical example:
```
&system
  seed='symm'
  beta=2000.d0,
  mu=2.0
  spectra_calc = .true.
/
&control
  use_lehman = .false.
  trace_only = .false.
  ff_only    = .false.
  npade = 40

  fast_calc  = .true.
  nnu=1
  emin = 0.0001
  emax = 0.0
/
```

Keywords are:

* seed : seedname of hr.dat file.

* beta : inverse temperature to be used.

* mu   : system chemical potential

* spectra_calc : If .true., the calculation is done for real frequencies, as specified by _nnu_, _emin_ and _emax_, otherwise, first _nnu_ Matsubara frequencies are calculated.

* nnu/emin/emax : For real-frequency calculations, a frequency mesh from _emin_ to _emax_ with equally-distributed _nnu_ frequencies are calculated; for Matsubara frequency calculations, the first _nnu_ Matsubara frequencies are calculated.

* use_lehman : If .true., the calculation is done using Lehman representation, otherwise it is done by summing up Green's functions

* trace_only : If .true., full matrix is not generated. Only compatible with use_lehman=.true.

* ff_only    : If .true., FC / CF / CC blocks are not calculated. Only compatible with trace_only=.false.

* npade      : Number of poles used in Pad\'{e} summation. Only compatible with use_lehman=.false.

* fast_calc  : If .true. use fast algorithm. More memory is required than fast_calc=.false.