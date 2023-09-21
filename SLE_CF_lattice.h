//
//  SLE_CF_lattice.h
//  LE-LP
//
//  Created by moulay abdella chkifa on 10/13/21.
//  Copyright © 2021 moulay abdella chkifa. All rights reserved.
//

#ifndef SLE_CF_lattice_h
#define SLE_CF_lattice_h

#include "SLE_alias.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif


/*
***************************************************
So-called Hadamard matrices used in the paper.
These are square matrices of dim 1,2,4,8,16,...
defined by recursion
***************************************************
 */
mat_d Hadamard_matrix (size_t n){
    int N = (int)pow(2,n);
    mat_d Hn(N,vec_d(N,1));
    if (n==0){return Hn;}
    else{
        mat_d X = Hadamard_matrix (n-1) ;
        for (int i=0; i<N/2; i++) {
            for (int j=0; j<N/2; j++) {
                Hn[i][j]     = X[i][j];
                Hn[i][j+N/2] = X[i][j];
                Hn[i+N/2][j] = X[i][j];
                Hn[i+N/2][j+N/2] = -X[i][j];
            }
        }
        return Hn;
    }
}
/*
**********************************************************
Hadamard matrices shrunk to have determinant 1/N,
then multiplied by 2 to be enumerated in [-1,1]^d
dimension of matrix is d=2^n and scaling factor is N
The matrix B is the inverse of transpose of A
*********************************************************
 */
std::vector<mat_d> Hadamard_lattice_matrices (size_t n, double N){
    int d = (int)pow(2,n);
    
    vec_d row (d);
    
    mat_d A (d,row);
    mat_d AT(d,row); // transpose of A
    mat_d B (d,row);
    mat_d BT(d,row); // transpose of B
    
    double shrinkFactor = 0.5*pow(N,1./d)*sqrt(d);//For hadamard matrices
    mat_d Hn = Hadamard_matrix(n);
    for (int i=0; i<d; i++) {
        for (int j=0; j<d; j++) {
            double Aij = Hn[i][j]/shrinkFactor;
            A[i][j] = Aij;
            AT[j][i] = Aij;
            double Bij = Hn[i][j]*shrinkFactor/d;
            B[i][j] = Bij ;
            BT[j][i] = Bij;
        }
    }
    return std::vector<mat_d>{A,AT,B,BT};
}

/*
 Square Tchebechev-Vandermonde Lattices of any kind.
 The general term is cos (j theta_i) for j=0,...,d-1
 and d angles theta_0,theta_1,...
 Output matrices are normalized for our convenience
 */
std::vector<mat_d> TV_lattice_matrices (size_t d,
                                        std::vector<double>& angles,
                                        double alphaA,
                                        double alphaB){
        
    std::vector<int> coefsA (d,2); coefsA[0]=1;
    std::vector<int> coefsB (d,1);
    
    vec_d row (d);
    mat_d A (d,row);
    mat_d AT(d,row); // transpose of A
    mat_d B (d,row);
    mat_d BT(d,row); // transpose of B
    vec_d k (d);
    
    for (int i=0; i<d; i++) {
        for (int j=0; j<d; j++) {
            double cos_ji = cos(j*angles[i]);
            if (fabs(cos_ji)<0.00000000001)
                cos_ji=0.0;
            double Aij = coefsA[j]*cos_ji/alphaA;
            A[i][j]  = Aij;
            AT[j][i] = Aij;
            double Bij = coefsB[j]*cos_ji*alphaB;
            B[i][j]  = Bij;
            BT[j][i] = Bij;
        }
    }
    return std::vector<mat_d>{A,AT,B,BT};
}
/*
 *******************************************************************
 Tchebechev Lattice of the orthogonal kind
 Tchebechev-Vandermonde matrices is shrunk to have determinant
 1/N, then multiplied by 2 to be enumerated in [-1,1]^d
 dimension of matrix is d and scaling factor is N
 The matrix B is the inverse of transpose of A
 *******************************************************************
 */
std::vector<mat_d> Tchebychev_lattice_matrices (size_t d, size_t N){
    std::vector<double> angles(d);
    for (int i=0; i<d; i++)
        angles[i] = (i+0.5)*M_PI/d;

    double shrinkFactor = 0.5*pow(N/sqrt(2),1./d) *sqrt(2.*d);
    return TV_lattice_matrices(d, angles, shrinkFactor, shrinkFactor/d);
}

#endif /* SLE_CF_lattice_h */
