program amgx_miniapp
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscmat.h>
!#include "solver_amgx.h"
   use petscksp
   use petscmat
   use iso_c_binding
   implicit none
   character*(128)  mat_file, rhs_file, which
   character(kind=c_char,len=128) solver_type
   PetscViewer      fd
   PetscBool flg
   PetscInt its
   PetscReal norm, solution_norm
   integer ierr, nprocs, color, key, rank
   integer amgxcomm
   KSP :: ksp
   Mat :: A
   Vec :: lhs, rhs, x
   type(c_ptr) :: amgx_struct
   integer(c_long_long) :: mataddr, rhs_addr, lhs_addr
   double precision :: t1, t2
   PetscScalar      none

    interface 
       subroutine allocate_amgx_struct(amgxstruct)  bind(C, name='allocate_amgx_struct') 
          use iso_c_binding
          type(c_ptr):: amgxstruct
       end subroutine allocate_amgx_struct
       subroutine init_amgx(amgxstruct, amgxcomm) bind(C, name='init_amgx')
           use iso_c_binding
           type(c_ptr) :: amgxstruct
           integer :: amgxcomm
       end subroutine init_amgx
       subroutine solve_amgx(amgxstruct, mataddr, rhs_addr, lhs_addr) bind(C, name='solve_amgx')
           use iso_c_binding
           type(c_ptr) :: amgxstruct
           integer(c_long_long) :: mataddr, rhs_addr, lhs_addr
       end subroutine solve_amgx
       subroutine getiters_amgx(amgxstruct, niter) bind(C, name='getiters_amgx')
           use iso_c_binding
           type(c_ptr) :: amgxstruct
           integer :: niter
       end subroutine getiters_amgx
       subroutine terminate_solver_amgx(amgxstruct) bind(C, name='terminate_solver_amgx')
           use iso_c_binding
           type(c_ptr) :: amgxstruct
       end subroutine terminate_solver_amgx

    end interface

    none = -1.0
    call PetscInitialize('petscrc',ierr)
    if (ierr .ne. 0) then
      print*,'Unable to initialize PETSc'
      stop
    endif

   !Read in matrix and RHS
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,&
      PETSC_NULL_CHARACTER,'-mat',mat_file,flg,ierr);CHKERRA(ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,&
      PETSC_NULL_CHARACTER,'-rhs',rhs_file,flg,ierr);CHKERRA(ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,&
      PETSC_NULL_CHARACTER,'-which',which,flg,ierr);CHKERRA(ierr)


    call MPI_Comm_Size(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr)


   


    call PetscViewerBinaryOpen(MPI_COMM_WORLD,mat_file,FILE_MODE_READ,fd,ierr);
    call MatCreate(MPI_COMM_WORLD,A,ierr)
    call MatSetFromOptions(A, ierr)
    call MatLoad(A,fd,ierr)
    call PetscViewerDestroy(fd,ierr)
    call PetscViewerBinaryOpen(MPI_COMM_WORLD,rhs_file,FILE_MODE_READ,fd,ierr);
    call VecCreate(MPI_COMM_WORLD,rhs,ierr)
    call VecSetFromOptions(rhs, ierr)
    call VecLoad(rhs,fd,ierr)
    call VecDuplicate(rhs,lhs,ierr)
    call VecSetFromOptions(lhs, ierr)
    call VecDuplicate(rhs,x,ierr)
    call VecSetFromOptions(x, ierr)
    call PetscViewerDestroy(fd,ierr)

    if(trim(which).eq.'amgx')then
        call allocate_amgx_struct(amgx_struct)
        call init_amgx(amgx_struct, MPI_COMM_WORLD)
        mataddr = A%v
        rhs_addr = rhs%v
        lhs_addr = lhs%v
        t1 = MPI_Wtime()
        call solve_amgx(amgx_struct, mataddr, rhs_addr, lhs_addr)
        t2 = MPI_Wtime()
        call getiters_amgx(amgx_struct, its)
        call terminate_solver_amgx(amgx_struct)
    elseif(trim(which).eq.'petsc')then
        call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
        call KSPSetOperators(ksp,A,A,ierr)
        call KSPSetFromOptions(ksp,ierr)
        t1 = MPI_Wtime()
        call KSPSolve(ksp,rhs,lhs,ierr)
        t2 = MPI_Wtime()
        call KSPGetIterationNumber(ksp,its,ierr)
    endif
    call MatMult(A,lhs,x,ierr)
    call VecAXPY(x,none,rhs,ierr)
    call VecNorm(x,NORM_2,norm,ierr)
    call VecNorm(lhs,NORM_2,solution_norm,ierr)
   
    if(rank.eq.0) write(6,100) solution_norm, norm,its, t2-t1
100 format('Solution norm ',e14.8, ' Residual norm',e11.4,' iterations ',i5, ' time',e11.2)

    call PetscFinalize(ierr)
    call MPI_Finalize(ierr)

end program amgx_miniapp
