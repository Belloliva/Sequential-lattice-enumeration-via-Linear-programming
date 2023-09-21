//
//  SLE_LP_B.h
//  LE-LP
//
//  Created by moulay abdella chkifa on 9/19/21.
//  Copyright © 2021 moulay abdella chkifa. All rights reserved.
//
#ifndef SLE_LP_B_h
#define SLE_LP_B_h

#include <numeric>

#include "SLE_static.h"
#include "SLE_simplex.h"

inline void make_rhs_positive(mat_d& A, vec_d& b){
    size_t rowSize = A.size();
    size_t colSize = A[0].size();

    for(size_t i=0; i<rowSize; i++){
        int sign = (b[i]>=0)? 1: -1;
        b[i] = sign*b[i];
        for(int j=0; j<colSize; j++)
            A[i][j] = sign * A[i][j];
    }
}

/*
 Initialize linear programs in standard form whose constaints correspond to MT w = c, w>=0 where
 MT = [+Bi.T   | - Bi.T  ]
      [------------------]
      [   I    |    I    ]
and c=(a,1) where matrices Bi consist in i leading columns of B.
Bi.T is the transpose of matrix Bi. In the paper, these are the
linear programs we use in computing the fonction beta and which
we refer to as dual LPs using matrix B, see equation 24.
*/
linear_programming::LP_solver usingB_dual_LP(size_t i,
                                             const vec_d& a,
                                             const mat_d& BT){
    size_t dim = a.size();
    size_t colSize = 2*dim;
    size_t rowSize = i+dim;
    /*----------------------------------------------------------------------------
     Compute the matrix MT, the rhs b in constraints and the objective vector c
     ----------------------------------------------------------------------------*/
    mat_d MT(rowSize, vec_d(colSize, 0)); // transpose of M as in the paper
    for(int k=0; k<i; k++){
        std::copy(BT[k].begin(), BT[k].end(), MT[k].begin());
        std::copy(BT[k].begin(), BT[k].end(), MT[k].begin()+dim);
        std::for_each(MT[k].begin()+dim, MT[k].end(),[]( double &x){ x*=-1;});
    }
    
    for(int k=0; k<dim; k++){
        MT[i+k][k] = 1.;
        MT[i+k][dim+k] = 1.;
    }
    
    vec_d b(rowSize,1.);
    for(int k=0; k<i; k++){b[k] = a[k];};
    
    vec_d c(colSize,0.);
    std::copy(BT[i].begin(), BT[i].end(), c.begin());
    std::copy(BT[i].begin(), BT[i].end(), c.begin()+dim);
    std::for_each(c.begin()+dim,c.end(),[]( double &x){ x*=-1;});
    
    // make sure the vector b has >= entries and make the necessary change in the matrix MT
    make_rhs_positive(MT, b);
    /*---------------------------------------------------------------------*/
    linear_programming::LP_data pb_inf ( std::move(MT), std::move(b), std::move(c));
    linear_programming::LP_solver solver_inf (std::move(pb_inf));
    return solver_inf;
}

/*
 Initialize linear programs in standard form whose constaints correspond to MT w = c, w>=0 where
 MT = [+Bi.T   |    0    ]
      [------------------]
      [   I    |    I    ]
and c=((a+Bi.T 1)/2,1) where matrices Bi consist in i leading columns of B.
Bi.T is the transpose of matrix (Bi). In the paper, these are the linear
programs we refer to as simplified dual LPs using matrix B, see equation 41.
*/
linear_programming::LP_solver usingB_dual_LP_simple(size_t i,
                                                    const vec_d& a,
                                                    const mat_d& BT){
    size_t dim = a.size();
    size_t colSize = 2*dim;
    size_t rowSize = i+dim;
    /*----------------------------------------------------------------------------
     Compute the matrix MT, the rhs b in constraints and the objective vector c
     ----------------------------------------------------------------------------*/
    mat_d MT(rowSize, vec_d(colSize, 0)); // transpose of M as in the paper
    for(int k=0; k<i; k++)
        std::copy(BT[k].begin(), BT[k].end(), MT[k].begin());
    
    for(int k=0; k<dim; k++){
        MT[i+k][k] = 1.;
        MT[i+k][dim+k] = 1.;
    }
    
    vec_d b(rowSize,1.);
    for(int k=0; k<i; k++){
        double delta = std::accumulate(BT[k].begin(), BT[k].end(), 0.0, std::plus<double>());
        b[k] = 0.5*(a[k] +delta) ;
    };
    
    vec_d c(colSize,0.);
    std::copy(BT[i].begin(), BT[i].end(), c.begin());
    std::copy(BT[i].begin(), BT[i].end(), c.begin()+dim);
    std::for_each(c.begin()+dim,c.end(),[]( double &x){ x*=-1;});
        
    // make sure the vector b has >= entries and make the necessary change in the matrix MT
    make_rhs_positive(MT, b);
    /*---------------------------------------------------------------------*/
    linear_programming::LP_data pb_inf ( std::move(MT), std::move(b), std::move(c));
    linear_programming::LP_solver solver_inf (std::move(pb_inf));
    return solver_inf;
}

/*
Solve minimimze and maximize <b,w> under constraint MT w = c, w>=0,
We solve linear programs based on the matrix B as explained in the paper,
see equation 24 in the paper. BT is the transpose matrix of B.
 */

inline std::pair<double,double> beta_dual_LP(const mat_d& BT,
                                             size_t i,
                                             const vec_d& a){
    if (i==0){
        double beta0 = Norm_l1(BT[0]);
        return std::make_pair (-beta0,+beta0);
    }
    linear_programming::LP_solver solver_inf = usingB_dual_LP(i, a, BT);
        
    solver_inf.bypass_sanitize();
    solver_inf.set_initial_basis();
    solver_inf.run_simplex_phase_I();
    
    if (solver_inf.is_unfeasible ){return std::make_pair(1,0);}

    // otherwise keep going, solver_inf is feasible and solver_sup
    // will also be feasible since it shares the same constraints.
    
    linear_programming::LP_solver solver_sup = solver_inf;
    vec_d& c_sup = solver_sup.LP.c ;
    std::for_each(c_sup.begin(),c_sup.end(),[]( double &x){ x*=-1;});
    
    solver_inf.run_simplex_phase_II();
    solver_sup.run_simplex_phase_II();
    
    double lower_beta = +solver_inf.getOptimum();
    double upper_beta = -solver_sup.getOptimum();
    
    return std::make_pair(lower_beta,upper_beta);
}

/*
 Solve minimimze and maximize <b,w> under constraint MT w = c, w>=0,
 We solve linear programs based on the matrix B as explained in the paper,
 see equation 24 which is simplified in equation 41 in the paper.
 BT is the transpose matrix of B.
 */
inline std::pair<double,double> beta_dual_LP_simple(const mat_d& BT,
                                                    size_t i,
                                                    const vec_d& a){
    if (i==0){
        double beta0 = Norm_l1(BT[0]);
        return std::make_pair (-beta0,+beta0);
    }

    linear_programming::LP_solver solver_inf = usingB_dual_LP_simple(i, a, BT);
    
    solver_inf.bypass_sanitize();
    solver_inf.set_initial_basis();
    solver_inf.run_simplex_phase_I();
    
    if (solver_inf.is_unfeasible ){return std::make_pair(1,0);}

    // otherwise keep going, solver_inf is feasible and solver_sup
    // will also be feasible since it shares the same constraints.
    
    linear_programming::LP_solver solver_sup = solver_inf;
    vec_d& c_sup = solver_sup.LP.c ;
    std::for_each(c_sup.begin(),c_sup.end(),[]( double &x){ x*=-1;});
    
    solver_inf.run_simplex_phase_II();
    solver_sup.run_simplex_phase_II();
    
    double lower_beta = +solver_inf.getOptimum();
    double upper_beta = -solver_sup.getOptimum();
    
    return std::make_pair(lower_beta,upper_beta);
}

/*
 ******************************************************************************
 The main recursive function in the enumeration. The action in the inner-most is
 either to store the node or to evaluate an integrand. For now we only increment
 the number of enumerated lattice points by 1
 ******************************************************************************
 */

void enumerate_usingB_dual_LP_recursion(int j,
                                        vec_d& k,
                                        const mat_d& A,
                                        const mat_d& BT){
    
    size_t d = A.size();
    if (j==d){
        count_seen+=1;
        vec_d x = dot_product(A, k, d);
        if (belong_to_hypercube(x)){count_enumerated+=1;}
    }
    
    else{
        std::pair <double,double> pr;
        //pr = beta_dual_LP(BT, j, k); // first option, vanilla
        pr = beta_dual_LP_simple(BT, j, k); // second option, faster
        
        double lower_j = pr.first  -0.1;
        double upper_j = pr.second +0.1;
        int min_kj = ceil(lower_j);
        int max_kj = floor(upper_j);
        if (max_kj < min_kj)
            count_dead_ends +=1;

        for (int kj =min_kj ; kj <= max_kj; kj++) {
            k[j] = kj;
            enumerate_usingB_dual_LP_recursion(j+1, k, A, BT);
        }
    }
}


void enumerate_usingB_dual_LP(const mat_d& A, const mat_d& BT){
    size_t dim = A.size();
    vec_d k(dim);
    enumerate_usingB_dual_LP_recursion(0, k, A, BT);
}


/* --------------------------------------------------------------------------------
 Using the predefined Linear programs
----------------------------------------------------------------------------------*/

/*
 Initialize linear programs in standard form whose constaints correspond to MT w = c, w>=0 where
 MT = [+Bi.T |     0     ]
      [------------------]
      [   I    |    I    ]
and c=((a+Bi.T 1)/2,1) where matrices Bi consist in i leading columns of B.
Bi.T is the transpose of matrix Bi. In the paper, these are the linear
programs we refer to as dual LPs using matrix B, see equation 41 and
which we use in computing the fonction beta
*/
inline std::vector<linear_programming::LP_solver> initial_LP_B(const mat_d& BT, int mode){
    
    size_t dim = BT.size();
    std::vector<linear_programming::LP_solver> vec_LP(dim );
    size_t num_cols = 2*dim;
    size_t num_rows = dim;
    
    vec_d v(num_cols,0);
    for (int i=0; i<num_rows; i++) {
        mat_d BBi(dim+i,v);
        
        for (int j=0; j<i; j++)
            std::copy(BT[j].begin(), BT[j].end(), BBi[j].begin());
        
        for (int j=0; j< dim ; j++){
            BBi[i+j][j]=1;
            BBi[i+j][dim+j]=1;
        }
        
        vec_d bbi(dim+i,1.);
        
        vec_d cci(2*dim);
        std::copy     (BT[i].begin(), BT[i].end(), cci.begin());
        std::copy     (BT[i].begin(), BT[i].end(), cci.begin()+dim);
        std::for_each(cci.begin()+dim, cci.end(), []( double &x){ x*=-1;});
        
        if (mode==SIMPLEX_MAXIMIZE)
            std::for_each(cci.begin(), cci.end(), []( double &x){ x*=-1;});

        linear_programming::LP_data pb( std::move(BBi), std::move(bbi), std::move(cci));
        // bbi has no role, it will be changed throught the algorithm
        // We end here, every time we will change bbi and solve using simplex algorithm.
        linear_programming::LP_solver LP_i(std::move(pb));
        LP_i.bypass_sanitize();// LP is consistent
        LP_i.set_initial_basis();
        //LP_i.run_simplex_phase_I();
        vec_LP [i] = std::move(LP_i);
    }
    return vec_LP;
}


inline std::pair<double ,double> beta_dual_LP_Predefined(int i,
                                               const vec_d& rhs,
                                               linear_programming::LP_solver& solver_inf){
    mat_d& A_inf = solver_inf.LP.A;
    vec_d& b_inf = solver_inf.LP.b;
    std::copy(rhs.begin(), rhs.end(), b_inf.begin());
    make_rhs_positive(A_inf, b_inf);
    solver_inf.run_simplex_phase_I();
    
    if (solver_inf.is_unfeasible ){return std::make_pair(1,0);}

    // otherwise keep going, solver_inf is feasible and solver_sup
    // will also be feasible since it shares the same constraints.

    linear_programming::LP_solver solver_sup = solver_inf;
    vec_d& c_sup = solver_sup.LP.c;
    std::for_each(c_sup.begin(), c_sup.end(), []( double &x){ x*=-1;});
    
    
    solver_inf.run_simplex_phase_II();
    solver_sup.run_simplex_phase_II();

    double lower_beta = + solver_inf.getOptimum();
    double upper_beta = - solver_sup.getOptimum();
    
    return std::make_pair( lower_beta , upper_beta);
}


/*
 ******************************************************************************
 The main recursive function in the enumeration
 The action in the inner-most is either to
 store the node or to evaluate an integrand
 ******************************************************************************
 */

void enumerate_usingB_dual_LP_recursion_Predefined(int j,
                                                   vec_d& k,
                                                   const mat_d& A,
                                                   const mat_d& BT,
                                                   const std::vector<vec_d>& sum_BT,
                                                   std::vector<linear_programming::LP_solver>& LPs_inf){
    
    size_t d = A.size();
    if (j==d){
        count_seen+=1;
         vec_d x = dot_product(A, k, d);
        if (belong_to_hypercube(x)){
            count_enumerated+=1;
            //std::cout <<"number enumerated lattice points = " << count_enumerated << std::endl;
        }
    }
    
    else{
        std::pair <double,double> pr;
            
        // we don't use a reference, so we work with a new program everytime
        linear_programming::LP_solver solver_inf =  LPs_inf[j];

        vec_d rhs(j);
        for (int i=0; i<j; i++)
            rhs[i] = 0.5*(k[i]+sum_BT[j][i]);
        pr = beta_dual_LP_Predefined(j, rhs, solver_inf);
        
        double lower_j = pr.first -0.1;
        double upper_j = pr.second+0.1;
        int min_kj = ceil(lower_j);
        int max_kj = floor(upper_j);
        if (max_kj < min_kj)
            count_dead_ends +=1;
    
        for (int kj =min_kj ; kj <= max_kj; kj++) {
            k[j] = kj;
            enumerate_usingB_dual_LP_recursion_Predefined(j+1, k, A, BT,sum_BT, LPs_inf);
        }
    }
}



void enumerate_usingB_dual_LP_Predefined(const mat_d& A, const mat_d& BT){
    
    size_t dim = A.size();
    vec_d k(dim);
    std::vector<linear_programming::LP_solver> LPs_inf = initial_LP_B(BT,SIMPLEX_MINIMIZE);
    //std::vector<linear_programming::LP_solver> LPs_sup = initial_LP_B(BT,SIMPLEX_MAXIMIZE);
    // The simplex code only solves minimization linear programs
    // The parameter mode = SIMPLEX_MINIMIZE or SIMPLEX_MAXIMIZE in the initial_LP function is used
    // in order to transform a maximization into minimization.
    std::vector<vec_d> sum_BT (dim);
    // The vector of all Bi.T times (1 ... 1) for i=0,...,d-1
    // Bi.T is the transpose of Bi the matrix of i leading columns of B
    for (int j=1; j<dim; j++){
        sum_BT[j] = vec_d(j);
        std::copy(sum_BT[j-1].begin(), sum_BT[j-1].end(), sum_BT[j].begin());
        sum_BT[j][j-1] = std::accumulate(BT[j-1].begin(), BT[j-1].end(), 0.0, std::plus<double>());
    }
    enumerate_usingB_dual_LP_recursion_Predefined(0, k, A, BT,sum_BT, LPs_inf);
}


#endif /* SLE_LP_B_h */
