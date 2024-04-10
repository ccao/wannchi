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

### Description of IBZKPT file

### Description of QPOINTS file

### Description of RPA.inp file

### Description of wannchi.inp
