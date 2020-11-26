program amgx_miniapp
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscmat.h>
!#include "solver_amgx.h"
   use petscksp
   use petscmat
   use iso_c_binding
   use amgx_mod
   implicit none
   character*(128)  mat0_file, rhs0_file, mat1_file, rhs1_file, which
   character(kind=c_char,len=128) solver_type
   PetscViewer      fd
   PetscBool flg
   PetscInt its
   PetscReal norm, solution_norm
   integer ierr, nprocs, color, key, rank
   integer amgxcomm
   KSP :: ksp
   Mat :: A0, A1
   Vec :: lhs0, rhs0, x0
   Vec :: lhs1, rhs1, x1
   type(c_ptr) :: amgx_struct
   integer(c_long_long) :: mataddr, rhs_addr, lhs_addr
   double precision :: t1, t2
   PetscScalar      none


    none = -1.0
    call PetscInitialize('petscrc',ierr)
    if (ierr .ne. 0) then
      print*,'Unable to initialize PETSc'
      stop
    endif

   !Read in matrix and RHS
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,&
      PETSC_NULL_CHARACTER,'-mat0',mat0_file,flg,ierr);CHKERRA(ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,&
      PETSC_NULL_CHARACTER,'-rhs0',rhs0_file,flg,ierr);CHKERRA(ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,&
      PETSC_NULL_CHARACTER,'-mat1',mat1_file,flg,ierr);CHKERRA(ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,&
      PETSC_NULL_CHARACTER,'-rhs1',rhs1_file,flg,ierr);CHKERRA(ierr)

    call PetscOptionsGetString(PETSC_NULL_OPTIONS,&
      PETSC_NULL_CHARACTER,'-which',which,flg,ierr);CHKERRA(ierr)


    call MPI_Comm_Size(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr)
   


    call PetscViewerBinaryOpen(MPI_COMM_WORLD,mat0_file,FILE_MODE_READ,fd,ierr);
    call MatCreate(MPI_COMM_WORLD,A0,ierr)
    call MatSetFromOptions(A0, ierr)
    call MatLoad(A0,fd,ierr)
    call PetscViewerDestroy(fd,ierr)

    call PetscViewerBinaryOpen(MPI_COMM_WORLD,rhs0_file,FILE_MODE_READ,fd,ierr);
    call VecCreate(MPI_COMM_WORLD,rhs0,ierr)
    call VecSetFromOptions(rhs0, ierr)
    call VecLoad(rhs0,fd,ierr)
    call PetscViewerDestroy(fd,ierr)

    call VecDuplicate(rhs0,lhs0,ierr)
    call VecSetFromOptions(lhs0, ierr)
    call VecDuplicate(rhs0,x0,ierr)
    call VecSetFromOptions(x0, ierr)

    call PetscViewerBinaryOpen(MPI_COMM_WORLD,mat1_file,FILE_MODE_READ,fd,ierr);
    call MatCreate(MPI_COMM_WORLD,A1,ierr)
    call MatSetFromOptions(A1, ierr)
    call MatLoad(A1,fd,ierr)
    call PetscViewerDestroy(fd,ierr)

    call PetscViewerBinaryOpen(MPI_COMM_WORLD,rhs1_file,FILE_MODE_READ,fd,ierr);
    call VecCreate(MPI_COMM_WORLD,rhs1,ierr)
    call VecSetFromOptions(rhs1, ierr)
    call VecLoad(rhs1,fd,ierr)
    call VecDuplicate(rhs1,lhs1,ierr)
    call VecSetFromOptions(lhs1, ierr)
    call VecDuplicate(rhs1,x1,ierr)
    call VecSetFromOptions(x1, ierr)

    call PetscViewerDestroy(fd,ierr)

    if(trim(which).eq.'amgx')then
        call allocate_amgx_struct(amgx_struct)
        call init_amgx(amgx_struct, MPI_COMM_WORLD, 1)
        mataddr = A0%v
        rhs_addr = rhs0%v
        lhs_addr = lhs0%v
        t1 = MPI_Wtime()
        call set_amgx_new(amgx_struct, mataddr, rhs_addr, lhs_addr, MPI_COMM_WORLD)
        t2 = MPI_Wtime()
        if(rank.eq.0 ) write(*,*) 'AMGX Setup time',  t2-t1
        t1 = MPI_Wtime()
        call solve_amgx_new(amgx_struct)
        t2 = MPI_Wtime()
        if(rank.eq.0 ) write(*,*) 'AMGX Solve time',  t2-t1

        call MatCopy(A1, A0, SAME_NONZERO_PATTERN, ierr)
        call VecCopy(rhs1, rhs0, ierr)
        t1 = MPI_Wtime()
        call update_amgx_new(amgx_struct)
        t2 = MPI_Wtime()
        if(rank.eq.0 ) write(*,*) 'AMGX Update time',  t2-t1
        t1 = MPI_Wtime()
        call solve_amgx_new(amgx_struct)
        t2 = MPI_Wtime()
        if(rank.eq.0 ) write(*,*) 'AMGX Solve time',  t2-t1


        call terminate_solver_amgx(amgx_struct)
    elseif(trim(which).eq.'petsc')then
        call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
        call KSPSetFromOptions(ksp,ierr)
        call KSPSetOperators(ksp,A0,A0,ierr)
        t1 = MPI_Wtime()
        call KSPSolve(ksp,rhs0,lhs0,ierr)
        t2 = MPI_Wtime()
        if(rank.eq.0 ) write(*,*) 'PETSC Solve',  t2-t1
        if(rank.eq.0) write(*,*) "#####################################"
        call KSPSetOperators(ksp,A1,A1,ierr)
        t1 = MPI_Wtime()
        call KSPSolve(ksp,rhs1,lhs1,ierr)
        t2 = MPI_Wtime()
        if(rank.eq.0 ) write(*,*) 'PETSC Solve',  t2-t1
        call KSPGetIterationNumber(ksp,its,ierr)
    endif
    call MatMult(A0,lhs0,x0,ierr)
    call VecAXPY(x0,none,rhs0,ierr)
    call VecNorm(x0,NORM_2,norm,ierr)
    call VecNorm(lhs0,NORM_2,solution_norm,ierr)
   
    if(rank.eq.0) write(6,100) solution_norm, norm,its, t2-t1
100 format('Solution norm ',e14.8, ' Residual norm',e11.4,' iterations ',i5, ' time',f12.5)

    call PetscFinalize(ierr)
    !call MPI_Finalize(ierr)

end program amgx_miniapp
