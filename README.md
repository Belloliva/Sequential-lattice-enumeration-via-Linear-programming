# SLE-LP
C++ code for the paper "Lattice enumeration via linear programming" 

- Chkifa, M.A. Lattice enumeration via linear programming. Numer. Math. (2023). https://doi.org/10.1007/s00211-023-01376-6

structure of the code

  + timer: a class for a timer, to produce runtimes
  + SLE_alias: aliases for a vector of integer, vector of doubles, matrix of doubles (vector of rows)
  + SLE_static: static int variables used to count number of enumerated lattice points
  + SLE_linalg: few simple linear algebra functions
  + SLE_simplex: an implemenation of the simplex tableau algorithm
  + SLE_CF_lattice: functions for producing the generating matrices for Hadamard lattices and orthogonal Chebyshev-Frolov lattices
  + SLE_LP_A : enumeration functions based on the generating matrix A
  + SLE_LP_B : enumeration functions based on the generating matrix B of dual lattice
    
