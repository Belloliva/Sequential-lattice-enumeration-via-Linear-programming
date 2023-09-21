//
//  SLE_python_print.h
//  LE-LP
//
//  Created by moulay abdella chkifa on 9/19/21.
//  Copyright © 2021 moulay abdella chkifa. All rights reserved.
//

#ifndef SLE_python_print_h
#define SLE_python_print_h

#include <vector>
#include <algorithm>    // std::minmax
#include <iomanip>
#include "SLE_alias.h"


template <class container>
void print_to_python (const container& v){
    std::cout << std::fixed;
    std::cout << std::setprecision(1);
    std::cout << "[";
    size_t i;
    for (i=0; i<v.size()-1; i++)
        std::cout << v[i]  << "\t,";
    std::cout << v[i] << "]" ;
}


template <class container>
void print_to_python_np (const std::string& s, const container& v){
    std::cout << s << "np.array(";
    print_to_python (v);
    std::cout << ")\n";
}

template <class Mat_container>
int print_to_python_np_ (const Mat_container & A){
    std::cout << "np.array([\n";
    size_t i;
    for (i=0; i<A.size()-1; i++){
        print_to_python (A[i]);
        std::cout << ",\n";
    }
    print_to_python (A[i]);
    std::cout << "])\n";
    return 0;
}

void print_to_python_np (const std::string& s, const mat_d & A){
    std::cout << s;
    print_to_python_np_ (A);
}

#endif /* SLE_python_print_h */
