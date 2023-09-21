# SLE-LP
C++ code for the paper "Lattice enumeration via linear programming" to appear in numerische mathematik

structure of the code

  + timer: a class for a timer, to produce runtimes
  + SLE_alias: aliases for a vector of integer, vector of doubles, matrix of doubles (vector of rows)
  + SLE_static: static int variables used to count number of enumerated lattice points
  + SLE_linalg: few simple linear algeba functions
  + SLE_simplex: an implemenation of the simplex tableau alrgorithm
  + SLE_CF_lattice: functions for producing the generating matrices for Hadamard lattices and orthogonal Chebyshev-Frolov lattices
  + SLE_LP_A : enumeration functions based on the generating matrix A
  + SLE_LP_B : enumeration functions based on the generating matrix B of dual lattice
    
