# solver_miniapp
This is a minimalistic miniapp that allow to solve a linear system Ax=b by using either AMGX or PETSc
Dependencies:
PETSC, AMGX, AmgXWrapper
Execution
mpiexec -np 1 ./miniapp_amgx.x -mat petsc_mat -rhs petsc_rhs -which 'petsc' or 'amgx'
You also need a petscrc file or a amgx.conf file 
