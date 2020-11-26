MODULE amgx_mod
    interface 
       subroutine allocate_amgx_struct(amgxstruct)  bind(C, name='allocate_amgx_struct') 
          use iso_c_binding
          type(c_ptr):: amgxstruct
       end subroutine allocate_amgx_struct
       subroutine init_amgx(amgxstruct, amgxcomm, solver) bind(C, name='init_amgx')
           use iso_c_binding
           type(c_ptr) :: amgxstruct
           integer :: amgxcomm
           integer, value :: solver
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
       subroutine set_amgx_new(amgxstruct, mataddr, rhsaddr, lhsaddr, amgxcomm) bind(C, name='set_amgx_new')
           use iso_c_binding
           type(c_ptr) :: amgxstruct
           integer :: amgxcomm
           integer(c_long_long) :: mataddr, rhsaddr, lhsaddr
       end subroutine set_amgx_new
      subroutine solve_amgx_new(amgxstruct) bind(C, name='solve_amgx_new')
           use iso_c_binding
           type(c_ptr) :: amgxstruct
      end subroutine solve_amgx_new
      subroutine update_amgx_new(amgxstruct) bind(C, name='update_amgx_new')
           use iso_c_binding
           type(c_ptr) :: amgxstruct
      end subroutine update_amgx_new


    end interface
END MODULE
