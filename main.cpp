//
//  main.cpp
//  LE-LP
//
//  Created by moulay abdella chkifa on 9/18/21.
//  Copyright © 2021 moulay abdella chkifa. All rights reserved.
//

#include <cmath>
#include <numeric>
#include <iostream>

#include "SLE_CF_lattice.h"
#include "SLE_LP_A.h"
#include "SLE_LP_B.h"

#include "timer.h"

// A square matrix is reprsented as vector of rows, see the alias
// mat_d in SLE_alias. In order to enumerate all k such that
// -1 <= A k <= 1, you need to provide the matrix A and it transpose
// AT and use the fast function enumerate_usingA_dual_LP.
// In order to experiment with enumeration using the matrix B
// generating the dual lattice, you need to provide the tranpose BT of
// B (which is the inverse of A) and use the functions
//enumerate_usingB_dual_LP or enumerate_usingB_dual_LP_Predefined

// The matrices defined below were used in order
// to produce the numerical results of the paper

int main(int argc, const char * argv[]) {
    
    std::vector<float> runtimes;
    std::vector<size_t> numSeenPoints;
    std::vector<size_t> numEnumeratedPoints;
    
    std::vector<size_t> vec_dims(20);
    std::iota(vec_dims.begin(), vec_dims.end(), 8);
    
    std::vector<size_t> vec_m {16};
    //std::iota(vec_m.begin(), vec_m.end(), 16);
    
    for (size_t d:vec_dims) {
        for (size_t m:vec_m) {
            double N = pow(2,m);
            std::cout << "---------------------------------------------------" << std::endl;
            std::cout <<"-** "<< "d="<< d <<", m=" << m <<", N=" << N <<" **-" << std::endl;
            
            //std::vector<mat_d> AABB = Hadamard_lattice_matrices ( 4, N);
            std::vector<mat_d> AABB = Tchebychev_lattice_matrices ( d, N);
                        
            mat_d A = AABB[0]; // matrix A generating the lattice
            mat_d AT= AABB[1]; // transpose of A
            mat_d B = AABB[2]; // matrix B=transpose(inverse(A))
            mat_d BT= AABB[3]; // transpose of B, which is the inverse of A
            
            // initialize the static variables for every enumeration
            count_dead_ends=0;
            count_seen=0;
            count_enumerated=0;
            count_succes=0;

            // set the timer
            timer t ("algo_timer1");
            
            enumerate_usingA_dual_LP(A, AT);
            
            //enumerate_usingB_dual_LP(A, BT);
            //enumerate_usingB_dual_LP_Predefined(A, BT);
            
            std::cout << "Number of encountred dead-ends: " << count_dead_ends << std::endl;
            std::cout << "Number of seen lattice points: "       << count_seen << std::endl;
            std::cout << "Number of enumerated lattice points: " << count_enumerated << std::endl;
            //std::cout << "Number of feasabe LPs: " << count_feasability << std::endl;

            runtimes.push_back(t.timer_value());
            numSeenPoints.push_back(count_seen);
            numEnumeratedPoints.push_back(count_enumerated);
        }
    }
    
    std::cout << "##################################################" <<std::endl;
    std::cout << std::endl;

    print_to_python_np("dimensions d= ", vec_dims);
    print_to_python_np("runtimes = ", runtimes);
    print_to_python_np("numSeenPoints = ", numSeenPoints);
    print_to_python_np("numEnumeratedPoints = ", numEnumeratedPoints);
    
    return 0;
    
}
