#include "petscksp.h"
#include "petscmat.h"
#include "petscsys.h"
#include "petscvec.h"
#include "petscviewer.h"

#include "AmgXSolver.hpp"
#include "petscToCSR.hpp"

// Run with: CUDA_VISIBLE_DEVICES=0 mpirun -n <NRANKS> ./miniapp-all-mps.x <solver> <mat_files> <vec_files>
// e.g. for AMGX - CUDA_VISIBLE_DEVICES=0 mpirun -n <NRANKS> ./miniapp-all-mps.x amgx matrices/mat_poisson* matrices/vec_poisson*
//      for PETSc -  CUDA_VISIBLE_DEVICES=0 mpirun -n <NRANKS> ./miniapp-all-mps.x petsc matrices/mat_poisson* matrices/vec_poisson*

#include <mpi.h>

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, "petscrc", NULL);
  PetscGetArgs(&argc, &argv);

  int rank, nranks, amgx_rank, amgx_nranks;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  MPI_Comm amgx_comm;
  MPI_Comm_split(MPI_COMM_WORLD, rank, 0, &amgx_comm);

  MPI_Comm_rank(amgx_comm, &amgx_rank);
  MPI_Comm_size(amgx_comm, &amgx_nranks);

  printf("rank %d of %d, rank %d of %d in AMGX\n", rank, nranks, amgx_rank,
         amgx_nranks);

  // Input format is mat_1 mat_2 mat_3 vec_1 vec_2 vec_3 etc.

  Mat A;
  Vec rhs, lhs, x;
  PetscViewer fdmat, fdvec;
  KSP ksp;

  MatCreate(amgx_comm, &A);
  MatSetFromOptions(A);

  VecCreate(amgx_comm, &lhs);
  VecCreate(amgx_comm, &rhs);
  VecCreate(amgx_comm, &x);
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

  int div = num_inputs / nranks;
  int rem = num_inputs % nranks;

  int beg_inputs = rank * div + std::min(rank, rem);
  int end_inputs = (rank + 1) * div + std::min(rank + 1, rem);

  printf("rank %d, beg = %d, end = %d\n", rank, beg_inputs, end_inputs);

  KSPCreate(amgx_comm, &ksp);
  KSPSetFromOptions(ksp);

  char* amgx_conf_file = getenv("AMGX_CONF_FILE");
  if (amgx_conf_file == NULL) amgx_conf_file = (char*)"amgx.conf";

  AmgXSolver solver;
  solver.initialize(amgx_comm, "dDDI", amgx_conf_file);

  std::vector<double> times;
  std::vector<double> norms;
  std::vector<double> sol_norms;
  std::vector<std::string> mat_names;
  std::vector<std::string> vec_names;
  std::vector<int> iters;

  bool setA = false;

  PetscSynchronizedPrintf(MPI_COMM_WORLD, "rank %d num inputs = %d\n", rank,
                          end_inputs - beg_inputs);
  PetscSynchronizedFlush(MPI_COMM_WORLD, stdout);

  for (int i = beg_inputs; i < end_inputs; ++i) {
    char** mat = inputs + i;
    char** vec = inputs + i + num_inputs;

    PetscViewerBinaryOpen(amgx_comm, *mat, FILE_MODE_READ, &fdmat);
    PetscViewerBinaryOpen(amgx_comm, *vec, FILE_MODE_READ, &fdvec);

    MatLoad(A, fdmat);
    VecLoad(rhs, fdvec);
    VecDuplicate(rhs, &lhs);
    VecDuplicate(rhs, &x);

    int its = 0;

    double t0 = MPI_Wtime();

    PetscInt nRowsLocal, nRowsGlobal, nNz;
    const PetscInt* colIndices = nullptr;
    const PetscInt* rowOffsets = nullptr;
    PetscScalar* values;
    double* amgx_lhs;
    double* amgx_rhs;

    petscToCSR(amgx_comm, A, nRowsLocal, nRowsGlobal, nNz, rowOffsets,
               colIndices, values, amgx_lhs, amgx_rhs, lhs, rhs);

    if (!setA) {
      solver.setA(nRowsGlobal, nRowsLocal, nNz, rowOffsets, colIndices, values,
                  nullptr);
      setA = true;
    }
    else {
      solver.updateA(nRowsLocal, nNz, values);
    }

    solver.solve(amgx_lhs, amgx_rhs, nRowsLocal);
    solver.getIters(its);

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

  PetscPrintf(MPI_COMM_WORLD, "%40s %40s %20s %20s %20s %20s %8s\n", "mat-file",
              "vec-file", "iters", "time", "solution-norm", "residual-norm",
              "rank");

  double avg_time = 0.0;

  for (int i = 0; i < (end_inputs - beg_inputs); ++i) {
    PetscSynchronizedPrintf(MPI_COMM_WORLD,
                            "%40s %40s %20d %20.5f %20.5e %20.5e %8d\n",
                            mat_names[i].c_str(), vec_names[i].c_str(),
                            iters[i], times[i], sol_norms[i], norms[i], rank);
    avg_time += times[i];
  }
  PetscSynchronizedFlush(MPI_COMM_WORLD, stdout);

  double all_avg_time = 0.0;
  MPI_Reduce(&avg_time, &all_avg_time, 1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  PetscPrintf(MPI_COMM_WORLD, "Average time across ranks = %f\n",
              all_avg_time / nranks / num_inputs);

  solver.finalize();

  PetscFinalize();
  return 0;
}
