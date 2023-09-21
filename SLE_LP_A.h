//
//  SLE_LP_A.h
//  LE-LP
//
//  Created by moulay abdella chkifa on 9/19/21.
//  Copyright © 2021 moulay abdella chkifa. All rights reserved.
//
#ifndef SLE_LP_A_h
#define SLE_LP_A_h

#include <numeric>

#include "SLE_static.h"
#include "SLE_simplex.h"

/*
###########################################################################
// Using the generating matrix A for enumeration
############################################################################
 */

/*
Ak is a vector that consists in A1k1,A2k2,...,Adkd where Aj is matrix of leading j
columns of A and kj is the section of leading j coordinates of k.
Every time we set a new coordinate k_i of the integer vector k, we need to reflects the
change in Ak, Aiki is equal Ajkj with j=i-1 plus k_i a_i where a_i is the ith column of A
*/
inline void Update_partial_Ak (size_t i, const vec_d& k, std::vector<vec_d>& Ak, const mat_d& AT){
    if (i==-1) return;
    const vec_d& Ai = AT[i];
    double ki = k[i];
    
    if (i==0)
        std::transform(Ai.begin(),
                       Ai.end(),
                       Ak[0].begin(),
                       [ki](double x){return ki*x;});
    else
        std::transform(Ai.begin(),
                       Ai.end(),
                       Ak[i-1].begin(),
                       Ak[i].begin(),
                       [ki](double x,double y){return ki*x+y;});
}

/*
Initialize linear programs in standard form whose constaints correspond to
M = [(Ai->).T -(Ai->).T] and c=e1 where matrices (Ai->) consist in (d-i)
trailing columns of A and (Ai->).T is the transpose of matrix (Ai->). In the
paper, these are the linear programs we refer to as dual LPs using matrix A,
see equation 34 of the paper. We also initialize the simplex tableaux using
phase I simplex, we already know that the constaints are well defined and
consistant, hence the use of bypass_sanitize
*/
inline std::vector<linear_programming::LP_solver> initial_LP(const mat_d& AT, int mode){
    
    size_t dim = AT.size();
    // we initialize all the tableaux needed during the execution of enumeration algorithm
    std::vector<linear_programming::LP_solver> vec_LP(dim );
    size_t num_cols = 2*dim;
    size_t num_rows = dim;
    
    vec_d v(num_cols,0);
    for (int i=0; i<num_rows; i++) {
        mat_d AAi(num_rows-i,v);// Matrix M transpose in equation 34 of the paper
        for (int j=i; j<num_rows; j++) {
            vec_d& curr_row = AAi[j-i];
            std::copy     (AT[j].begin(), AT[j].end(), curr_row.begin());
            std::transform(AT[j].begin(), AT[j].end(), curr_row.begin()+dim, [](double  x) {return -x;});
            // if we are initializing tableaux for maximization, we multiply by -1 the current row
            if (mode==SIMPLEX_MAXIMIZE)
                std::for_each(curr_row.begin(), curr_row.end(), []( double &x){ x*=-1;});
        }
        vec_d bbi(num_rows-i,0.);
        bbi[0]=1;
        
        linear_programming::LP_data pb(std::move(AAi), std::move(bbi), vec_d(2*dim,0));
        // Since we only execute phase I simplex, the objective vector is irrelevant at this stage.
        
        linear_programming::LP_solver LP_i(std::move(pb));
        LP_i.bypass_sanitize();
        LP_i.set_initial_basis();
        LP_i.run_simplex_phase_I();
        
        vec_LP [i] = std::move(LP_i);
    }
    return vec_LP;
}

/*
Intsantiate the linear programs derived in equation 34 of the paper. This
function is only called by alpha_dual_LP
*/
inline std::pair<linear_programming::LP_solver,linear_programming::LP_solver>
usingA_dual_LP(size_t i, const vec_d& a, const mat_d& A, const mat_d& AT){
    size_t dim = a.size();
    size_t colSize = 2*dim;
    size_t rowSize = dim-i;
    
    
    mat_d MT_inf(rowSize, vec_d(colSize, 0)); // transpose of M as in the paper
    for(int k=0; k<rowSize; k++){
        std::copy(AT[i+k].begin(), AT[i+k].end(), MT_inf[k].begin());
        std::copy(AT[i+k].begin(), AT[i+k].end(), MT_inf[k].begin()+dim);
        std::for_each(MT_inf[k].begin()+dim,MT_inf[k].end(),[]( double &x){ x*=-1;});
    }
    
    mat_d MT_sup(rowSize, vec_d(colSize, 0)); // transpose of M as in the paper
    for(int k=0; k<rowSize; k++){
        std::copy(AT[i+k].begin(), AT[i+k].end(), MT_sup[k].begin());
        std::copy(AT[i+k].begin(), AT[i+k].end(), MT_sup[k].begin()+dim);
        std::for_each(MT_sup[k].begin(),MT_sup[k].begin()+dim,[]( double &x){ x*=-1;});
    }
    
    
    //should initialize the arrays b and c here
    vec_d b_inf(rowSize,0.); b_inf[0]=1;
    vec_d c_inf(colSize,1.);
    vec_d Aa = dot_product(A, a, i); //0 if i=0
    for (int k=0; k<dim; k++) { c_inf [k] -=Aa[k]; c_inf [dim+k] +=Aa[k];}
    
    
    vec_d b_sup{b_inf};
    vec_d c_sup{c_inf};
        
    linear_programming::LP_data pb_inf(std::move(MT_inf), std::move(b_inf), std::move(c_inf));
    linear_programming::LP_solver solver_inf (std::move(pb_inf));
    
    linear_programming::LP_data pb_sup(std::move(MT_sup), std::move(b_sup), std::move(c_sup));
    linear_programming::LP_solver solver_sup (std::move(pb_sup));
    
    
    return std::make_pair(solver_sup, solver_inf);
    // solver_sup is used to derive lower alpha
    // solver_inf is used to derive upper alpha
}

/*
Solve the dual LPs maximaze{M.T w = -c, w>=0} - b w and minimize{M.T w = c, w>=0}  b w
derived at equation 34 of the paper.
A is the generating matrix, and AT is the transpose matrix of A.
Here we intsantiate new linear programs everytime we need to solve them,
which make the code slow and inefficient. This function is here for experiementing
it is not used.
*/
inline std::pair<double,double> alpha_dual_LP(const mat_d& A , const mat_d& AT, size_t i, const vec_d& a){
    
    auto pr = usingA_dual_LP(i, a, A, AT);
    linear_programming::LP_solver solver_sup = pr.first;
    linear_programming::LP_solver solver_inf = pr.second;
    
    solver_sup.bypass_sanitize();
    solver_sup.set_initial_basis();
    solver_sup.run_simplex_phase_I();
    solver_sup.run_simplex_phase_II();
    
    solver_inf.bypass_sanitize();
    solver_inf.set_initial_basis();
    solver_inf.run_simplex_phase_I();
    solver_inf.run_simplex_phase_II();

        
    double lower_alpha = -solver_sup.getOptimum();
    double upper_alpha = +solver_inf.getOptimum();
    return std::make_pair( lower_alpha, upper_alpha);
}

/*
Using the chained simlpex algorithm, the idea is simple, we solve LPs using the simplex
tableau algorithm, the algorithm converge and we leave the simplex tableaux at their final
states, so they can be used for warm starting the algorthim on a subsequent request.
*/
inline std::pair<double ,double> alpha_chained_dual_LP(int i,
                                                       const std::vector<vec_d>& vec_Aa,
                                                       linear_programming::LP_solver& solver_sup,
                                                       linear_programming::LP_solver& solver_inf){
    
    size_t dim = vec_Aa.size();
    vec_d Aa(dim,0);
    if (i>=1){ Aa = vec_Aa[i-1];}
    // The objective vector is (1-Aa | 1+Aa) for minimisation (and
    // maximization if turned into minimization) linear programs
    // derived in equation 34 of the paper, see vector b in equation 32.

    vec_d& c_sup = solver_sup.LP.c;
    vec_d& c_inf = solver_inf.LP.c;

    std::for_each(c_sup.begin(), c_sup.end(), []( double &x){ x=+1;} );
    std::transform(c_sup.begin(),     c_sup.begin()+dim, Aa.begin(), c_sup.begin()    , std::minus<double>());
    std::transform(c_sup.begin()+dim, c_sup.end()      , Aa.begin(), c_sup.begin()+dim, std::plus<double>());
    
    std::copy(c_sup.begin(), c_sup.end(), c_inf.begin());

    // Now we solve linear program which we know are feasible and the simplex
    // tableau algorithm was already initialized for their constaints, so no
    // need to initialize a basis or run simplex phase I

    solver_sup.optimum=0;
    solver_sup.run_simplex_phase_II();
    if (solver_sup.is_unbounded)// in this case no need to solve solver_inf it will also be unbouded
        return std::make_pair( 1 , 0);
    
    solver_inf.optimum=0;
    solver_inf.run_simplex_phase_II();
    
    double lower_alpha = - solver_sup.getOptimum();
    double upper_alpha = + solver_inf.getOptimum();
        
    return std::make_pair( lower_alpha , upper_alpha);
}


void enumerate_usingA_dual_LP_recusion(int j,
                                       vec_d& k,
                                       std::vector<vec_d>& Ak,
                                       const mat_d& A,
                                       const mat_d& AT,
                                       std::vector<linear_programming::LP_solver>& LPs_sup,
                                       std::vector<linear_programming::LP_solver>& LPs_inf){
    
    size_t d = A.size();
    if (j==d){
        count_seen+=1;
        Update_partial_Ak (d-1, k, Ak, AT); // so that the last element in Ak is A times k
        vec_d& x=Ak[d-1];//can also define x as follow//vec_d x = dot_product(A, k, d);
        if (belong_to_hypercube(x))
            count_enumerated+=1;
    }
    else{
        Update_partial_Ak (j-1, k, Ak, AT); // update the vector of products Ak up to j
        // If we use references, we implement the FLE-LP algorithm, otherwise
        // we are simply using the plain algorithm LE-LP-A with mutualized simplex
        linear_programming::LP_solver& solver_sup =  LPs_sup[j];
        linear_programming::LP_solver& solver_inf =  LPs_inf[j];
                
        /* Then we compute the [min_kj,min_kj] where kj should belong */
        // If we use alpha_dual_LP, we are simply implementing LE-LP-A
        // If we use alpha_chained_dual_LP, we make profit of accelaration based on
        // the idea of the chained simplex algorithm as explained in the paper
        
        // auto pr = alpha_dual_LP(A, AT, j, k);
        auto pr = alpha_chained_dual_LP( j, Ak, solver_sup, solver_inf);

        double lower_j = pr.first ;
        double upper_j = pr.second ;
        int min_kj, max_kj;
        // min_kj = ceil(lower_j); max_kj = floor(upper_j);
        // The previous two lines of code were commented and replaced by the next
        // lines as guard against numerical precision issues, like when upper_j
        // is 2.000001 but the simplex solver gives back 1.999999 hence min_kj is 1
        // but must be 2. We think the value 0.1 is adequat and safe, although in the
        // paper we suggest the safest value of 1
        
        if (upper_j==my_INFINITY){min_kj=1;max_kj=0;}
        // unbounded LP, so the next loop is not even entered and a dead end is counted
        else{
            double corrParar = 0.01; // a correction parameter against numerical errors
            min_kj = ceil(lower_j - corrParar);
            max_kj = floor(upper_j+ corrParar);
        }
        if (max_kj < min_kj) {count_dead_ends++;}

        for (int kj =min_kj ; kj <= max_kj; kj++) {
            k[j] = kj;
            enumerate_usingA_dual_LP_recusion(j+1, k, Ak, A, AT, LPs_sup, LPs_inf);
        }
    }
}

void enumerate_usingA_dual_LP(const mat_d& A, const mat_d& AT){
    
    size_t dim = A.size();
    vec_d k(dim);
    std::vector<vec_d> Ak(dim,k);
    //these are A1k1,A2k2,...,Adkd where Aj is matrix of leading j
    // columns of A and kj is the section of j leading coordinates of k
    // we first initialize all the tableau needed for simplex method
    std::vector<linear_programming::LP_solver> LPs_sup = initial_LP(AT,SIMPLEX_MAXIMIZE);
    std::vector<linear_programming::LP_solver> LPs_inf = initial_LP(AT,SIMPLEX_MINIMIZE);
    // The simplex code only solves minimization linear programs
    // The parameter mode = SIMPLEX_MAXIMIZE or SIMPLEX_MINIMIZE in the initial_LP
    // function is used in order to transform a maximization into minimization.
    enumerate_usingA_dual_LP_recusion(0, k, Ak, A, AT, LPs_sup, LPs_inf);
}

#endif /* SLE_LP_A_h */
