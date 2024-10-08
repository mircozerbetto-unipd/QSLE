### Contents

1.  [Important Information](#1-important-information)
2.  [Citation Details](#2-citation-details)
3.  [Distribution License](#3-distribution-license)
4.  [Prerequisites](#4-prerequisites)
5.  [Installation](#5-installation)
6.  [Software guide](#6-software-guide)

---

### 1. Important information

QSLE-v1.0 is the first program to simulate a system using the Quantum-Stochastic Liouville Equation.
For the instruction reguarding how to use the software, please read the documentation in
the doc folder.

This program is based on:
1) Armadillo: C++ Library for Linear Algebra & Scientific Computing  
Copyright 2008-2024 Conrad Sanderson (https://conradsanderson.id.au)  
Copyright 2008-2016 National ICT Australia (NICTA)  
Copyright 2017-2024 Data61 / CSIRO  
2) spline.h : a simple cubic spline interpolation library without external dependencies
Copyright (C) 2011, 2014, 2016, 2021 Tino Kluge (ttk448 at gmail.com)

Authors:
  * Riccardo Cortivo - riccardo.cortivo@phd.unipd.it
  * Mirco Zerbetto   - mirco.zerbetto@unipd.it

---

### 2: Citation Details

Please cite the following papers if you use QSLE software in your research.  
Citations are useful for the continued development and maintenance of the library.

  * Jonathan Campeggio, Riccardo Cortivo and Mirco Zerbetto.
    A multiscale approach to coupled nuclear and electronic dynamics. I. Quantum-stochastic Liouville equation in natural internal coordinates
    The Journal of Chemical Physics, Vol. 158, pp. 244104, 2023.

  * Riccardo Cortivo, Jonathan Campeggio and Mirco Zerbetto.
    A multiscale approach to coupled nuclear and electronic dynamics. II. Exact and approximated evaluation of nonradiative transition rates
    The Journal of Chemical Physics, Vol. 158, pp. 244105, 2023.
    
  * Jonathan Campeggio, Antonino Polimeno and Mirco Zerbetto.
    DiTe2: Calculating the diffusion tensor for flexible molecules
    The Journal of Computational Chemistry, Vol. 40, pp. 697, 2018.

  * Conrad Sanderson and Ryan Curtin.  
    Armadillo: a template-based C++ library for linear algebra.  
    Journal of Open Source Software, Vol. 1, No. 2, pp. 26, 2016.  
  
  * Conrad Sanderson and Ryan Curtin.  
    Practical Sparse Matrices in C++ with Hybrid Storage and Template-Based Expression Optimisation.  
    Mathematical and Computational Applications, Vol. 24, No. 3, 2019.

---

### 3: Distribution License

This program is free software; you can redistribute it and/or modify it under the 
terms of the GNU General Public License as published by the Free Software Foundation;
either version 2 of the License, or any later version. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this 
program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street,
Fifth Floor, Boston, MA  02110-1301, USA.

Armadillo is licensed under the Apache License, Version 2.0 (the "License").
A copy of the License is included in the "LICENSE.txt" file (in the QSLE-v1.0/lib/armadillo-12.8.0 folder).

---

### 4: Prerequisites

QSLE program can be installed only on Linux systems via Cmake, with or without root access.

The CMake tool can be downloaded from https://www.cmake.org 
or (preferably) installed using the package manager on your system.

C++11 is required (it is is fully supported by GCC 4.8+ and Clang 3.4+).

OpenMP 3.1+ is required, even if you only want to use the serial version of the program.

The following libraries have to be available as shared libraries before building the QSLE executable:
- OpenBLAS
- LAPACK
- ARPACK
- SuperLU : only 5.2.x, 5.3.x, or 6.0.x (the last one must be compiled with default integer size (32 bits))
It is also necessary to install the corresponding development files for each library.
For example, when installing the `libopenblas` package, also install the `libopenblas-dev` package.

---

### 5: Installation

To install QSLE-v1.0, just unpack the tgz package file

    tar zxf QSLE-v1.0.tgz
    
step into the QSLE-v1.0 folder and run the build.sh script

    cd QSLE-v1.0
    ./build.sh

then type 1 and press enter. The information about the compiling procedure will
be printed out on the terminal.

If you want to clean up all the the built packages and the executable, run again the build.sh script

    ./build.sh

then type 2 and press enter.

---

### 6: Software guide

For more details reguarding the installation procedure and to understand how to use QSLE-v1.0,
please read the QSLE_User_Guide.pdf manual located in the doc folder. 
In addition, you can find other information in the articles suggested in the Section 2 and in
an upcoming article that describes this software in details.




