//
//  SLE_linalg.h
//  LE-LP
//
//  Created by moulay abdella chkifa on 9/15/23.
//  Copyright © 2023 moulay abdella chkifa. All rights reserved.
//

#ifndef linalg_h
#define linalg_h

#include "SLE_alias.h"

// dot product of first j coordinates of two vectors
double  dot_product (const vec_d& v1,
                     const vec_d& v2,
                     size_t j){
    double ans= 0 ;
    // if i==0, we don't even execute the loop,
    // which agrees with convention of the paper.
    for (int k=0; k<j; k++)
        ans += v1[k]*v2[k];
    return ans;
}

// matrix-vector product of first j column of a matrix M
// with first j coordinates of the vector v
// a matrix is represented as a vector of its rows
vec_d dot_product (const mat_d& M,
                   const vec_d& v,
                   size_t j){
    size_t d = M.size();
    vec_d w(d);
    for (int i=0; i<d; i++)
        w[i] = dot_product (M[i], v, j);
    return w;
}

// matrix-matrix product of two squares matrices A and B
// a matrix is represented as a vector of its rows
// The matrix BT transpose of B is provided instead of B
// in order to perform fast computation using the previously
// defined dot_product function

mat_d mat_product(const mat_d& A, const mat_d& BT){
    size_t d = A.size();
    mat_d C(d,vec_d(d));
    for (int i=0; i<d; i++) {
        for (int j=0; j<d; j++) {
            C[i][j] = dot_product(A[i], BT[j], d);
        }
    }
    return C;
}


double Norm_l1(const vec_d& v){
    double ans = 0.;
    for (int j=0; j<v.size(); j++)
        ans += fabs(v[j]);
    return ans;
}

double Norm_l1(const vec_d& v, size_t i){
    double ans = 0.;
    for (int j=0; j<i; j++)
        ans += fabs(v[j]);
    return ans;
}

double Norm_l2(const vec_d& v){
    double ans = 0.;
    for (int j=0; j<v.size(); j++)
        ans += v[j]*v[j];
    return sqrt(ans);
}

double Norm_l2(const vec_d& v, size_t i){
    double ans = 0.;
    for (int j=0; j<i; j++)
        ans += v[j]*v[j];
    return sqrt(ans);
}


// For a vector to belong to [-1,1]^d
bool belong_to_hypercube(const vec_d& v){
    return std::all_of(v.begin(), v.end(), [](double x){return (abs(x)<=1);});
}

#endif /* linalg_h */
