static char help[] = "Miniapp";
static char petscrc[] = "petscrc";
#include "stdio.h"
#include "AmgXSolver.hpp"
#include <petscmat.h>
#include "petscksp.h"
#include <mpi.h>
#include "petscToCSR.hpp"

struct amgx_struct{
    AmgXSolver solver;
    const PetscInt      *colIndices=nullptr,
                        *rowOffsets=nullptr;

    PetscScalar         *values;
    PetscInt            nRowsLocal,
                        nRowsGlobal,
                        nNz;

    double              *lhs,
                        *rhs; 

};



extern "C" void terminate_solver_amgx(amgx_struct *& amgx){
    amgx->solver.finalize();
    return;
}
extern "C" void init_amgx(amgx_struct *& amgx,  MPI_Fint *amgxcomm, int which){
    MPI_Comm Ccomm = MPI_Comm_f2c(*amgxcomm);	
    printf("INIT AMGX %d %d\n", which);
    if(which==1)
       amgx->solver.initialize(Ccomm,"dDDI","amgx.poisson");
    if(which==2)
       amgx->solver.initialize(Ccomm,"dDDI","amgx.ampere");

    return;
}
extern "C" void set_amgx(amgx_struct *& amgx, Mat  *A){
    double t1, t2, t3, t4;
    PetscReal norm;
    t1 = MPI_Wtime();
    //amgx->solver.setA(amgx->A);
    amgx->solver.setA(*A);
    t2 = MPI_Wtime();
    //VecNorm(*lhs,  NORM_2, &norm);
    return;
}
extern "C" void set_amgx_new(amgx_struct *& amgx, Mat  *A,Vec *rhs_petsc, Vec  *lhs_petsc, MPI_Fint *amgxcomm){
    MPI_Comm Ccomm = MPI_Comm_f2c(*amgxcomm);	
    petscToCSR(Ccomm, *A, amgx->nRowsLocal, amgx->nRowsGlobal, amgx->nNz, amgx->rowOffsets, amgx->colIndices, amgx->values, amgx->lhs, amgx->rhs, *lhs_petsc, *rhs_petsc);
    amgx->solver.setA(amgx->nRowsGlobal, amgx->nRowsLocal, amgx->nNz, amgx->rowOffsets, amgx->colIndices, amgx->values, nullptr);
    //petscToCSR(MPI_COMM_WORLD, *A, amgx->nRowsLocal, amgx->nRowsGlobal, amgx->nNz, amgx->rowOffsets, amgx->colIndices, amgx->values, amgx->lhs, amgx->rhs, *lhs_petsc, *rhs_petsc);
    //amgx->solver.setA(amgx->nRowsGlobal, amgx->nRowsLocal, amgx->nNz, amgx->rowOffsets, amgx->colIndices, amgx->values, nullptr);
    //amgx->solver.solve(amgx->lhs, amgx->rhs, amgx->nRowsLocal);
    //printf("AAAAAAA %d %d %d\n", amgx->nRowsLocal, amgx->nRowsGlobal, amgx->nNz);

}
extern "C" void update_amgx_new(amgx_struct *& amgx){
    amgx->solver.updateA(amgx->nRowsLocal, amgx->nNz, amgx->values);
    return;

}
extern "C" void solve_amgx_new(amgx_struct *& amgx){
    amgx->solver.solve(amgx->lhs, amgx->rhs, amgx->nRowsLocal);
    return;
}
extern "C" void solve_amgx(amgx_struct *& amgx, Vec *rhs, Vec *lhs){
    double t1, t2, t3, t4;
    PetscReal norm;
    t3 = MPI_Wtime();
    amgx->solver.solve(*lhs, *rhs);
    t4 = MPI_Wtime();
    VecNorm(*lhs,  NORM_2, &norm);
    return;
}
extern "C" void allocate_amgx_struct(amgx_struct *& amgx){
   amgx = new amgx_struct();
   return;
}
/*
extern "C" void solve_amgx(amgx_struct *& amgx, Mat  *A, Vec *rhs, Vec *lhs){
    double t1, t2, t3, t4;
    PetscReal norm;
    t1 = MPI_Wtime();
    //amgx->solver.setA(amgx->A);
    amgx->solver.setA(*A);
    t2 = MPI_Wtime();
    t3 = MPI_Wtime();
    amgx->solver.solve(*lhs, *rhs);
    t4 = MPI_Wtime();
    printf("Setup_time Solve_time %f %f\n",t2-t1,t4-t3);
    VecNorm(*lhs,  NORM_2, &norm);
    return;
}
*/

   /*
extern "C" void  prepare_for_solver(amgx_struct *& amgx, char *input_mat, char *input_rhs){
    PetscViewer fd;
//amgx = new amgx_struct();
   //amgx_struct *amg = allocate_amgx_struct(MPI_COMM_WORLD); 
    MPI_Comm_size(amgx->Ccomm, &(amgx->size));
    MPI_Comm_rank(amgx->Ccomm, &(amgx->rank));
    PetscViewerBinaryOpen(amgx->Ccomm,input_mat,FILE_MODE_READ,&fd);
    MatCreate(amgx->Ccomm,&(amgx->A));
    MatSetFromOptions(amgx->A);
    MatLoad(amgx->A,fd);
    PetscViewerDestroy(&fd);
    PetscViewerBinaryOpen(amgx->Ccomm,input_rhs,FILE_MODE_READ,&fd);
    VecCreate(amgx->Ccomm,&(amgx->rhs));
    VecSetFromOptions(amgx->rhs);
    VecLoad(amgx->rhs,fd);
    PetscViewerDestroy(&fd);
    VecSetFromOptions(amgx->rhs);
    VecDuplicate(amgx->rhs, &(amgx->lhs));
    VecSetFromOptions(amgx->lhs);
    return ;
}

extern "C" void solver_amgx_(char *input_mat, char *input_rhs, MPI_Fint *amgxcomm){
    double t1, t2, t3, t4;
    int rank, size;
    amgx_struct *amgx = new amgx_struct();
    
    Mat A;
    Vec rhs, lhs;
    AmgXSolver solver;
    PetscReal norm;
    PetscViewer fd;
    MPI_Comm Ccomm;
    amgx->Ccomm = MPI_Comm_f2c(*amgxcomm);
    MPI_Comm_rank(amgx->Ccomm, &rank);
    MPI_Comm_size(amgx->Ccomm, &size);
    printf("BBBBB %d %d\n", rank, size);
    PetscViewerBinaryOpen(amgx->Ccomm,input_mat,FILE_MODE_READ,&fd);
    MatCreate(amgx->Ccomm,&(amgx->A));
    MatSetFromOptions(amgx->A);
    MatLoad(amgx->A,fd);
    PetscViewerDestroy(&fd);
    PetscViewerBinaryOpen(amgx->Ccomm,input_rhs,FILE_MODE_READ,&fd);
    VecCreate(amgx->Ccomm,&rhs);
    VecSetFromOptions(rhs);
    VecLoad(rhs,fd);
    PetscViewerDestroy(&fd);
    VecSetFromOptions(rhs);
    VecDuplicate(rhs, &lhs);
    VecSetFromOptions(lhs);

    amgx->solver.initialize(amgx->Ccomm,"dDDI","amgx.info");
    t1 = MPI_Wtime();
    amgx->solver.setA(amgx->A);
    amgx->solver.solve(lhs, rhs);
    t2 = MPI_Wtime();
    VecNorm(lhs,NORM_2,&norm);
    printf("AMGx norm %f %f\n", norm, t2-t1);
    amgx->solver.finalize();
    return;
}
else{
       t1 = MPI_Wtime();
       PetscViewerBinaryOpen(amgxcomm,ampere_mat,FILE_MODE_READ,&fd);
       MatCreate(amgxcomm,&A);
       MatSetFromOptions(A);
       MatLoad(A,fd);
       PetscViewerDestroy(&fd);
       PetscViewerBinaryOpen(amgxcomm,ampere_rhs,FILE_MODE_READ,&fd);
       VecCreate(amgxcomm,&rhs);
       VecSetFromOptions(rhs);
       VecLoad(rhs,fd);
       PetscViewerDestroy(&fd);
       VecSetFromOptions(rhs);
       VecDuplicate(rhs, &lhs);
       VecSetFromOptions(lhs);
       VecScatterCreateToZero(rhs, &ctx, &rhs_local);

       KSPCreate(amgxcomm, &ksp);
       KSPSetOperators(ksp,A,A);
       KSPSetFromOptions(ksp);
       KSPSolve(ksp,rhs,lhs);

       VecNorm(lhs,NORM_2,&norm);
       t2 = MPI_Wtime();
       printf("PETSc_norm  time %f %f\n", norm, t2-t1);
   
       //VecDuplicate(rhs, &lhs);
       //VecSetFromOptions(lhs);

   }
   MPI_Barrier(MPI_COMM_WORLD);
   t4 = MPI_Wtime();
   printf("Total time %f \n", t4-t3);
   PetscFinalize();
   return 0;
}
*/
