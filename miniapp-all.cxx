#include "petscksp.h"
#include "petscmat.h"
#include "petscsys.h"
#include "petscvec.h"
#include "petscviewer.h"

#include "AmgXSolver.hpp"
#include "petscToCSR.hpp"

// Run with: ./miniapp-all.x <solver> <mat_files> <vec_files>
// e.g. for AMGX - ./miniapp-all.x amgx matrices/mat_poisson* matrices/vec_poisson*
//      for PETSc -  ./miniapp-all.x petsc matrices/mat_poisson* matrices/vec_poisson*

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, "petscrc", NULL);
  PetscGetArgs(&argc, &argv);

  // Input format is mat_1 mat_2 mat_3 vec_1 vec_2 vec_3 etc.

  Mat A;
  Vec rhs, lhs, x;
  PetscViewer fdmat, fdvec;
  KSP ksp;

  MatCreate(MPI_COMM_WORLD, &A);
  MatSetFromOptions(A);

  VecCreate(MPI_COMM_WORLD, &lhs);
  VecCreate(MPI_COMM_WORLD, &rhs);
  VecCreate(MPI_COMM_WORLD, &x);
  VecSetFromOptions(lhs);
  VecSetFromOptions(rhs);
  VecSetFromOptions(x);

  std::string solver_type(argv[1]);
  argc--;
  argv++;

  if ((argc - 1) % 2) {
    fprintf(stderr, "error, unbalanced number of inputs\n");
    PetscFinalize();
    return 1;
  }

  char** inputs = argv + 1;
  int num_inputs = (argc - 1) / 2;

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetFromOptions(ksp);

  char* amgx_conf_file = getenv("AMGX_CONF_FILE");
  if (amgx_conf_file == NULL) amgx_conf_file = (char*)"amgx.conf";

  AmgXSolver solver;
  solver.initialize(MPI_COMM_WORLD, "dDDI", amgx_conf_file);

  std::vector<double> times;
  std::vector<double> norms;
  std::vector<double> sol_norms;
  std::vector<std::string> mat_names;
  std::vector<std::string> vec_names;
  std::vector<int> iters;

  for (int i = 0; i < num_inputs; ++i) {
    char** mat = inputs + i;
    char** vec = inputs + i + num_inputs;

    PetscViewerBinaryOpen(MPI_COMM_WORLD, *mat, FILE_MODE_READ, &fdmat);
    PetscViewerBinaryOpen(MPI_COMM_WORLD, *vec, FILE_MODE_READ, &fdvec);

    MatLoad(A, fdmat);
    VecLoad(rhs, fdvec);
    VecDuplicate(rhs, &lhs);
    VecDuplicate(rhs, &x);

    int its = 0;

    double t0 = MPI_Wtime();

    if (solver_type == "petsc") {
      KSPSetOperators(ksp, A, A);
      KSPSolve(ksp, rhs, lhs);
      KSPGetIterationNumber(ksp, &its);
    }
    else if (solver_type == "amgx") {
      PetscInt nRowsLocal, nRowsGlobal, nNz;
      const PetscInt* colIndices = nullptr;
      const PetscInt* rowOffsets = nullptr;
      PetscScalar* values;
      double* amgx_lhs;
      double* amgx_rhs;

      petscToCSR(MPI_COMM_WORLD, A, nRowsLocal, nRowsGlobal, nNz, rowOffsets,
                 colIndices, values, amgx_lhs, amgx_rhs, lhs, rhs);

      if (i == 0) {
        solver.setA(nRowsGlobal, nRowsLocal, nNz, rowOffsets, colIndices,
                    values, nullptr);
      }
      else {
        solver.updateA(nRowsLocal, nNz, values);
      }

      solver.solve(amgx_lhs, amgx_rhs, nRowsLocal);
      solver.getIters(its);
    }
    else {
      fprintf(stderr, "solver %s not recognised\n", solver_type.c_str());
    }

    double t1 = MPI_Wtime();

    double norm, solution_norm;

    MatMult(A, lhs, x);
    VecAXPY(x, -1.0, rhs);
    VecNorm(x, NORM_2, &norm);
    VecNorm(lhs, NORM_2, &solution_norm);

    times.push_back(t1 - t0);
    norms.push_back(norm);
    sol_norms.push_back(solution_norm);
    mat_names.push_back(std::string(*mat));
    vec_names.push_back(std::string(*vec));
    iters.push_back(its);

    PetscViewerDestroy(&fdmat);
    PetscViewerDestroy(&fdvec);
  }

  PetscPrintf(MPI_COMM_WORLD, "%40s %40s %20s %20s %20s %20s\n", "mat-file",
              "vec-file", "iters", "time", "solution-norm", "residual-norm");

  for (int i = 0; i < num_inputs; ++i) {
    PetscPrintf(MPI_COMM_WORLD, "%40s %40s %20d %20.5f %20.5e %20.5e\n",
                mat_names[i].c_str(), vec_names[i].c_str(), iters[i], times[i],
                sol_norms[i], norms[i]);
  }

  solver.finalize();

  PetscFinalize();
  return 0;
}
