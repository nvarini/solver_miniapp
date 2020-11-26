/*
 * Copyright (c) 2020, NVIDIA CORPORATION. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 * \file petscToCSR.hpp
 * \brief Prototypes of functions in petscToCSR.cpp
 * \author Matt Martineau (NVIDIA mmartineau@nvidia.com)
 * \date 2020-08-05
 */

#pragma once



// PETSc
# include <petscmat.h>
# include <petscvec.h>



PetscErrorCode petscToCSR(
    MPI_Comm comm,
    Mat &A,
    PetscInt& nRowsLocal,
    PetscInt& nRowsGlobal,
    PetscInt& nNz,
    const PetscInt*& rowOffsets,
    const PetscInt*& colIndices,
    PetscScalar*& values,
    double*& lhs,
    double*& rhs,
    Vec& lhs_petsc,
    Vec& rhs_petsc);
