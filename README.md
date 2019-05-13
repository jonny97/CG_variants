# CG_variants
This is the repo that contains the files for Jiannan Jiang's course project for both EE227C and CS267.

# CS 267 project files:
   stored in folder cori_setup, which does not contain the installation of poski and petsc, which it relies on.
   the scipy runs are in the notebook called "Python_implmentation_CACG_and_various_testing"
   ## test data:
      All the testing matrix can be found at https://sparse.tamu.edu/
   ## How to use
      To use the code, you have to install pOSKI (https://bebop.cs.berkeley.edu/poski/) first.
      To use multiple processors, you have to configure your machine with Petsc first, so the code can use MPI.     
   ## Tuned parameters for project report:
      Section 5: all the testing uses OneD Data Partitioning, two threads.
   
# EE227C project files:
   This projects demostrates a "smoothing" variant of CG that performs better than original CG in harder problems.
   



