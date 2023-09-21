//
//  SLE_simplex.h
//  LE-LP
//
//  Created by moulay abdella chkifa on 9/18/21.
//  Copyright © 2021 moulay abdella chkifa. All rights reserved.
//

#ifndef SLE_simplex_h
#define SLE_simplex_h

#include <vector>
#include <algorithm>    // std::minmax
#include "SLE_alias.h"
#include "SLE_linalg.h"
#include "SLE_python_print.h"

#ifndef my_INFINITY
#define my_INFINITY (1000000000)
#endif

static int print_simplex_num_iter = 0;

namespace linear_programming {

#define EPS 1E-8
#define DEBUG_SIMPLEX 0
#define PREPROCESS_SIMPLEX 0

/*
 Self explanatory functions, not used yet
 */

inline bool is_zero     (double x){ return fabs(x) <= EPS;}
inline bool is_not_zero (double x){ return fabs(x)  > EPS;};


inline bool is_zero     (vec_d& x){
    return std::all_of(x.begin(), x.end(), [](double x){return is_zero(x);});
}

inline bool is_not_zero (vec_d& x){
    return std::any_of(x.begin(), x.end(), [](double x){return is_not_zero(x);});
}

inline bool is_not_greater_zero (vec_d& x){
    return std::any_of(x.begin(), x.end(), [](double x){return x < -EPS;});
}

inline bool is_not_smaller_zero (vec_d& x){
    return std::any_of(x.begin(), x.end(), [](double x){return x > EPS;});
}

/*
 *************************************************
 Pivot a matrix given a pivot row and column.
 *************************************************
 */
inline void standard_pivoting (mat_d & A,
                               size_t row_index,
                               size_t col_index){
    
    vec_d& pivot_row = A[row_index];
    double p = pivot_row[col_index];
    
    if (fabs (p) <= EPS){
        std::cout <<"can not pivot, the pivot entry in too small"<< std::endl;
        exit(1);
    }
    
    /* normlize pivot row with value of pivot: */
    std::for_each(pivot_row.begin(), pivot_row.end(), [p](double & x){x/=p;});
    
    /* Do for all rows other than pivot row: */
    int i=0;
    for (vec_d& curr_row:A) {
        double p2 = curr_row[col_index];
        if ( i != row_index && fabs (p2) > EPS ) {
            std::transform(curr_row.begin(),
                           curr_row.end(),
                           pivot_row.begin(),
                           curr_row.begin(),
                           [p2](double x,double y){return x- p2*y ;});
        }
        i++;
    }
}
/*
 ******************************************************
 Pivot a matrix and a rhs given a pivot row and column
 *******************************************************
 */
inline void standard_pivoting_with_rhs (mat_d & A,
                                        vec_d & rhs,
                                        size_t row_index,
                                        size_t col_index){
    
    vec_d& pivot_row = A[row_index];
    double p = pivot_row[col_index];
    
    if (fabs (p) <= EPS){
        std::cout <<"can not pivot, the pivot entry in too small"<< std::endl;
        exit(1);
    }
    
    /* normlize pivot row with value of pivot: */
    std::for_each(pivot_row.begin(), pivot_row.end(), [p](double & x){x/=p;});
    rhs[row_index] /= p;
    
    /* Do for all rows other than pivot row: */
    int i=0;
    for (vec_d& curr_row:A) {
        double p2 = curr_row[col_index];
        if ( i != row_index and fabs (p2) > EPS ) {
            std::transform(curr_row.begin(),
                           curr_row.end(),
                           pivot_row.begin(),
                           curr_row.begin(),
                           [p2](double& x,double& y){return x- p2*y ;});
            rhs [i] -= p2 * rhs [row_index];
        }
        i++;
    }
}

/*
 **************************************************************
 Pivot a matrix given a pivot row and column.
 The pivoting is only performed on rows passed the pivot row
 **************************************************************
 */
inline void lower_pivoting (mat_d & A,
                            size_t row_index,
                            size_t col_index){
    size_t m = A.size();
    
    vec_d& pivot_row = A[row_index];
    double p = pivot_row[col_index];
    
    if (fabs (p) <= EPS){
        std::cout <<"can not pivot, the pivot entry in too small"<< std::endl;
        exit(1);
    }
    
    std::for_each(pivot_row.begin(), pivot_row.end(), [p](double & x){x/=p;});
    
    for (size_t i=row_index+1; i<m; i++) {
        vec_d& curr_row = A[i] ;
        double p2 = curr_row[col_index];
        
        if (fabs (p2) > EPS ){
            std::transform(curr_row.begin(),
                           curr_row.end(),
                           pivot_row.begin(),
                           curr_row.begin(),
                           [p2](double x,double y){return x- p2*y ;});
        }
    }
}

/*
 **************************************************************
 Pivot a matrix and a rhs given a pivot row and column.
 The pivoting is only performed on rows passed the pivot row
 **************************************************************
 */
inline void lower_pivoting_with_rhs (mat_d & A,
                                     vec_d & rhs,
                                     size_t row_index,
                                     size_t col_index){
    size_t m = A.size();
    
    vec_d& pivot_row = A[row_index];
    double p = pivot_row[col_index];
    
    if (fabs (p) <= EPS){
        std::cout <<"can not pivot, the pivot entry in too small"<< std::endl;
        exit(1);
    }
    
    std::for_each(pivot_row.begin(), pivot_row.end(), [p](double & x){x/=p;});
    rhs[row_index] /= p;
    
    for (size_t i=row_index+1; i<m; i++) {
        vec_d& curr_row = A[i] ;
        double p2 = curr_row[col_index];
        
        if (fabs (p2) > EPS ){
            std::transform(curr_row.begin(),
                           curr_row.end(),
                           pivot_row.begin(),
                           curr_row.begin(),
                           [p2](double x,double y){return x- p2*y ;});
            rhs [i] -= p2 * rhs[row_index];
        }
    }
}


/*
 ******************************************************************
 reduce a matrix to row echelon form and return its rank
 ******************************************************************
 */
size_t reduce_to_row_echelon_form(mat_d& A) {
    size_t num_rows = A.size();
    size_t num_cols = A[0].size();
    size_t rank = 0;
    
    size_t h = 0; /* Initialization of the pivot row */
    size_t k = 0; /* Initialization of the pivot column */
    
    while (h < num_rows and k < num_cols){
        /* Find the k-th pivot: */
        size_t i_max = h;
        double max_val = fabs(A[i_max][k]);
        size_t i = i_max+1;
        for (; i < num_rows; i++) {
            double curr_val = fabs(A[i][k]);
            if ( curr_val > max_val){i_max=i; max_val=curr_val;}
        }
        if (max_val < EPS) /* No pivot in this column, pass to next column */
            k ++;
        else{
            std::iter_swap(A.begin()+h, A.begin()+i_max);
            lower_pivoting(A, h, k);
            /* Increase pivot row and column */
            h++ ; k++ ; rank++;
        }
    }
    return rank;
}

/*
 ******************************************************************
 reduce a matrix to row echelon form and return its rank.
 An rhs is also reduced using the same operations
 ******************************************************************
 */
size_t reduce_to_row_echelon_form_with_rhs(mat_d& A, vec_d& rhs ) {
    size_t num_rows = A.size();
    size_t num_cols = A[0].size();
    size_t rank = 0;
    
    size_t h = 0; /* Initialization of the pivot row */
    size_t k = 0; /* Initialization of the pivot column */
    
    while (h < num_rows and k < num_cols){
        /* Find the k-th pivot: */
        size_t i_max = h;
        double max_val = fabs(A[i_max][k]);
        size_t i = i_max+1;
        for (; i < num_rows; i++) {
            double curr_val = fabs(A[i][k]);
            if ( curr_val > max_val){i_max=i; max_val=curr_val;}
        }
        if (max_val < EPS) /* No pivot in this column, pass to next column */
            k ++;
        else{
            std::iter_swap(A.begin()+h, A.begin()+i_max);
            std::iter_swap(rhs.begin()+h, rhs.begin()+i_max);
            lower_pivoting_with_rhs(A, rhs, h, k);
            /* Increase pivot row, pivot column and rank */
            h++ ; k++ ; rank++;
        }
    }
    return rank;
}

/*
 Class for the notion of basis in simplex algorithm, not used to its fullest yet
*/

class LP_basis{
    size_t num_rows;
    size_t num_cols;
    
    std::vector <size_t> basis_complement; // complement of the basis;
    std::vector <bool> rows_in; // flags for rows associated with basis
    
    std::vector <std::pair<size_t,size_t>> basis_rc; // bijection rows to columns;
    std::vector <std::pair<size_t,size_t>> basis_cr; // bijection columns to rows;
public:
    std::vector <int> basis; // indices of columns making the simplex basis;
    
public:
    LP_basis()=default;
    LP_basis(const LP_basis& other) = default;
    //~LP_basis()=default;
    
    LP_basis (size_t num_r,size_t num_c):basis{},num_rows{num_r}, num_cols{num_c},basis_rc{},basis_cr{}{
        basis = std::vector <int> (num_r,-1);
        //  rows_in = std::vector <bool> rows_in(num_rows,false);
        //  basis_complement= std::vector <size_t>(num_cols);
        //  std::iota(basis_complement.begin(), basis_complement.end(), 0);
    }
    
    bool contains (size_t j) const{
        return (std::find (basis.begin(),basis.end(),j) != basis.end());
    }
    
    size_t size() const {
        return std::count_if(basis.begin(),basis.end(), [](int x){return x>=0;});
    }
    
    
    // LP_basis (size_t dim):dimension{dim}, basis{std::vector <int>(dim,-1)}{}
    int& operator[](size_t j){ return basis[j];  }
    const int operator[](size_t j)const{ return basis[j];}
    
};

/*
Class for the data of a linear program in standard form
minimize c.x subject to constaints Ax=b, x>=0
*/

class LP_data {
public:
    mat_d  A;            // constraint matrix
    vec_d  b;            // right hand side
    vec_d  c;            // objective vector
    
public:
    LP_data()=default;
    //~LP_data()=default;
    
    LP_data(size_t rowSize,size_t colSize){
        A = mat_d(rowSize,vec_d(colSize));
        b = vec_d(rowSize);
        c = vec_d(colSize);
    }
    
    LP_data(const mat_d& A_,
            const vec_d& b_,
            const vec_d& c_):A{A_},b{b_},c{c_}{}
    
    LP_data(mat_d&& A_,
            vec_d&& b_,
            vec_d&& c_):A(std::move(A_)),b(std::move(b_)),c(std::move(c_)){}
    
    LP_data(const LP_data& other):A{other.A},b{other.b},c{other.c}{}// = default;
    LP_data( LP_data&& other)noexcept = default;
    
    LP_data& operator=(const LP_data &rhs)  = default;
    LP_data& operator=(LP_data &&rhs) noexcept = default;
    
    friend class LP_solver;
};


/*
Class for the linear program which is solved using revised simplex algorithm
The solver is designed for minimization
minimize c.x under the constraint Ax=b and x>=0
 */

class LP_solver{
public:
    LP_data LP;
    LP_basis basis;
    vec_d basic_solution;
    vec_d solution;               // unknowns
    double optimum;                  // optimum to be found
    bool is_well_defined    = false; // dimensions are consistent
    bool is_inconsistent    = false; // constraints are inconsistent
    bool is_unbounded       = false; // if well_defined and consistent but unbounded
    bool is_unfeasible      = false; // if well_defined, consistent but infeasibe
    bool is_solved          = false; // if well_defined, consistent, feasible and bounded, then it can be solved
    size_t number_redundancies = 0;
    
public:
    LP_solver()=default;
    //~LP_solver()=default;
    
    /*
     LP_solver(size_t rowSize,size_t colSize):basis{rowSize,colSize}{
        LP.A = mat_d(rowSize,vec_d(colSize));
        LP.b = vec_d(rowSize);
        LP.c = vec_d(colSize);
        solution = vec_d(colSize);
     }
     */
    
    LP_solver(const LP_data& LP_):LP{LP_},basis{LP.b.size(),LP.c.size()}{
        solution = vec_d (LP.c.size());
        //basis.basis = std::vector <int> (LP.b.size(),-1);
        //basis.basis = LP_basis (LP.b.size(),LP.c.size());
    }
    
    LP_solver( LP_data&& LP_):LP(std::move(LP_)),basis{LP.b.size(),LP.c.size()}{
        solution = vec_d (LP.c.size());
        //basis.basis  = LP_basis (LP.b.size(),LP.c.size());
    }
    
    
    
    // LP_solver(const LP_solver& other) = default;
    // LP_solver(const LP_solver& other) = default;
    // LP_solver& operator=(LP_solver &&rhs) noexcept = default;
    
    
    void check_sanity();
    void preprocess();
    void sanitize();
    void bypass_sanitize(){ // if we know the LP is well_defined, consistent, feasible and bounded
        is_well_defined    = true;  // dimensions are consistent
        is_inconsistent    = false; // constraint are consistent
        is_solved          = false; // not yet solved
    }
    
    void set_initial_basis();
    void run_simplex_phase_I();
    void run_simplex_phase_II();
    void run_simplex_phase(){ return; }
    
    vec_d getSolution(){
        if (is_solved) return solution;
        else{
            std::cout <<"the linear problem has not been solved yet"<< std::endl;
            exit(1);
        }
    }
    
    double getOptimum(){
        if (is_unbounded) return - my_INFINITY; // since we minimize
        if (is_solved) return - optimum; // we stick with the convention of revised simplex tableau
        else{
            std::cout <<"the linear problem has not been solved yet"<< std::endl;
            exit(1);
        }
    }
};
// Preprocessing:
// In order to check sanity and preprocess linear programs, we define helper
// functions which only act on matrix A and vectors b and c. Then we can call
// these functions for the linear program datas.
inline int cvn_check_sanity (mat_d& A, vec_d& b, vec_d& c);
inline int cvn_preprocess   (mat_d& A, vec_d& b, vec_d& c, vec_d &x);

//   First Make sure the linear program is consistent dimension-wise
inline int cvn_check_sanity(mat_d & A, vec_d & b, vec_d & c){
    size_t n = A[0].size();       // # of variables
    size_t m = A.size();          // # of inequalities
    
    if (!m or m != b.size() ){
        std::cout << "Wrong inputs: matrix A and vector b must have the same number of rows m>=1! \n";
        exit(1);
    }
    
    if (!n or n!= c.size()){
        std::cout << "Wrong inputs: matrix A and vector c must have the same number of columns n>=1! !\n";
        exit(1);
    }
    // then all rows must have the same dimension
    bool all_n = std::all_of (A.begin(), A.end(), [n](vec_d& row){ return row.size()==n;});
    if (!all_n){
        std::cout << "Wrong inputs: matrix rows must have the same dimension n>=1!\n";
        exit(1);
    }
    return 0;
}

// Preprocessing of Ax=b for use for revised simplex tableau:
//   1) Make sure that the r.h.s. (b) is non-negative.
//   2) Remove redundant constraints detected by linear dependency.
//   3) Solve the system if it is determined or over-determined, and check if it is consistent.
// Returns:
//   -1 : if the system Ax=b is inconsistent, or consistent has a unique solution x but x>=0 not fulfilled.
//   -2 : if the system is determined or over-determined and has unique
//        non-negative solution, which is stored in x.
//  >=0 : otherwise indicate the # of redundancies removed.

inline int cvn_preprocess(mat_d &A, vec_d &b, vec_d &c, vec_d &x){
    
    cvn_check_sanity(A,b,c);
    size_t m = A.size ();                 // # of constraints
    size_t n = A[0].size ();              // # of variables
    
    size_t rank = reduce_to_row_echelon_form_with_rhs(A, b);
    bool zero_b = std::all_of(b.begin()+rank, b.end(), [](const double& x){ return is_zero(x);});
    if(!zero_b) return -1; // inconsistent constraints
    
    if (m > n) {    // determined or overdetermined system (m>=n , so that rank <= n)
        // In this case, we solve a triagular linear system with A, b, x.
        // A is already in upper form, we use back-ward substitution
        std::for_each(x.begin(), x.end(), [](double& a){ a=0;});// make sure the entries of x are reset
        for (size_t i=rank-1; i > 0; i--) {
            vec_d& curr_row = A[i];
            x[i]= b[i] ;
            for (size_t j=i+1; j < rank; j++) {
                x[i] -= curr_row[j]*b[j];
            }
            if (x[i] < 0) return -1;// we exit if we encounter an entry that is negative
        }
        return -2;
    }
    int numRedundancies = (int)(m-rank);
    return numRedundancies;
}

void LP_solver::check_sanity(){
    int ss = cvn_check_sanity(LP.A, LP.b, LP.c);
    is_well_defined = (ss==0);
}

void LP_solver::preprocess(){
    int status = cvn_preprocess(LP.A, LP.b, LP.c, solution);
    if (status==-1)
        is_inconsistent=true;
    else if (status==-2)
        is_solved= true;
    else //  then prepare for simplex
        number_redundancies = status;
}

// By sanitize, we mean removes redundancies from well-defined, consistants constraints
void LP_solver::sanitize(){
    check_sanity();
    if(is_well_defined){
        preprocess();
        if(is_inconsistent){
            std::cout <<"the linear problem is inconsistent"<< std::endl;
            exit(1);
        }
        if(is_solved){
            std::cout <<"the linear problem is already solved"<< std::endl;
            exit(0);
        }
        // else, we delete redendencies and prepare the matrix/vector for simplex
        mat_d & A = LP.A; vec_d & b = LP.b;
        A.resize(A.size() - number_redundancies );
        b.resize(b.size() - number_redundancies );
        number_redundancies = 0;
    }
}

/*##########################################################
 // Solving using the simplex:
 // Self explanatory helper functions that are used for
 // phase I simplex algorithm
 ###########################################################*/
inline int is_unit_column (const mat_d & A, size_t c);
inline void concatenate_zeros(vec_d& v, size_t n);
inline void concatenate_kronecker(vec_d& v, size_t n, size_t j);

// Check if column c of A is a unit vector,
// meaning zero except for a 1 somewhere.
// Returns: -1 if not; non-zero row id otherwise.
inline int is_unit_column (const mat_d & A, size_t j) {
    int count  = 0;
    int row=0;
    double l1 = 0.;
    
    for (int r=0; r<A.size(); r++){
        l1+= fabs(A[r][j]);
        if (fabs(A[r][j] - 1.) < EPS){
            count++;
            row = r;
        }
    }
    return (count==1 and fabs(l1 - 1.) < EPS) ? row : -1;
}

inline void concatenate_zeros(vec_d& v, size_t n){
    v.resize(v.size() + n, 0.);
}

inline void concatenate_kronecker(vec_d& v, size_t n, size_t j){
    size_t old_size = v.size();
    v.resize(old_size + n, 0.);
    v[old_size + j] = 1.;
}

/*##########################################################
 // Solving using the simplex:
 // Self explanatory helper functions that are used used by
 // the simplex algorithm
 ##########################################################*/
inline int leaving_variable(mat_d & A, vec_d & b, size_t ev);
inline int entering_variable_simple_rule(const vec_d & c);
inline int entering_variable_dantzig_rule(const vec_d & c);
void canonicalize (mat_d & A, vec_d & b, vec_d & c, std::vector <int> & BasicVarR, double & obj);
bool pivoting     (mat_d & A, vec_d & b, vec_d & c, std::vector <int> & BasicVarR, double & obj, int max_iter);


// Find the index of a leaving variable
// in the simplex algorithm, given Ax=b
// and the index of the entering variable
inline int leaving_variable(mat_d & A,
                            vec_d & b,
                            size_t ev){
    // Find the leaving variable using ratio test.
    double minRatio=9999999999999.;
    int lv = -1, r=0;    // leaving variable, id'ed by row
    for (const vec_d& curr_row:A) {
        double piv = curr_row[ev];
        if (piv > EPS) {
            double ratio = b[r]/piv;
            if ( lv < 0 or ratio < minRatio ) {
                lv = r;
                minRatio = ratio;
            }
        }
        r++;
    }
    return lv;
    //if (lvr < 0) return true;      // unbounded
}

inline int entering_variable_simple_rule(const vec_d & c){
    // Find the entering variable (use the first one encountered).
    int ev = 0;                // id of the entering variable
    
    for (ev=0; ev<c.size(); ev++)
        if (c[ev] < -EPS) break;
    
    return (ev == c.size())? -1: ev;
}


inline int entering_variable_dantzig_rule(const vec_d & c){
    vec_d::difference_type ev =  std::min_element(c.begin(),c.end()) - c.begin();
    return (c[ev] >= -EPS )? -1: (int)ev;
}

// Convert a standard form of LP into canonical form:
// 1) The entries in the objective vector corresponding to basic variables are 0's.
// 2) The submatrix corresponding to basic variables is an identity matrix.
inline void canonicalize (mat_d & A,
                          vec_d & B,
                          vec_d & c,
                          std::vector <int> & BasicVarR,   // basic variable of each row
                          double & obj){              // objective value
    
    size_t m = A.size(), n = c.size();
    
    for (int r=0; r<m; r++) {
        int j = BasicVarR[r];          // col. that the basic variable is in
        
        
        if ( fabs(A[r][j] - 1.0) > EPS) {
            double p = A[r][j];
            for (int k=0; k<n; k++)
                A[r][k] /= p;
            B[r] /= p;
        }
        
        
        if (fabs(c[j]) > EPS) {
            double p = c[j];
            for (int k=0; k<n; k++)
                c[k] -= A[r][k] * p;
            obj -= B[r] * p;
        }
    }
}

// Pivoting iterations.
// Returns true if the problem is unbounded, false otherwise.
// This is the main function in the simplex algorithm
inline bool pivoting (mat_d & A,
                      vec_d & B,
                      vec_d & C,
                      std::vector <int> & BasicVarR,   // basic variable of each row
                      double & obj,
                      int max_iter){                   // objective value
    
    size_t m = A.size();
    size_t n = C.size();
    
    size_t count = 0;
    
    while (1 and count < max_iter){
        count++;
        
        int entering_var , leaving_var ;
        // Find the entering variable (use the one with minimum value of C).
        //entering_var = entering_variable_dantzig_rule(C);
        entering_var = entering_variable_simple_rule (C);
        if (entering_var<0 ){
            if( print_simplex_num_iter==1)
                std::cout << "simplex num iteration= " << count << std::endl;
            break;
        } // simplex converged
        
        leaving_var = leaving_variable( A, B, entering_var );
        if (leaving_var<0) return true; // problem unbounded
        
        // Leaving & entering variables identified, do pivoting.
        BasicVarR[leaving_var] = entering_var;
        standard_pivoting_with_rhs(A, B, leaving_var, entering_var);
        
        double p3 = C[entering_var];
        if ( fabs (p3) > EPS ) {
            vec_d& pivot_row =  A[leaving_var];
            
            for (int j=0; j<n; j++){
                C[j] -= p3*pivot_row[j];
            }
            /*std::transform(C.begin(),C.end(),pivot_row.begin(),C.begin(),[p3](double x,double y){return x- p3*y ;});*/
            obj -= p3 * B[leaving_var];
        }
        
        if (DEBUG_SIMPLEX) {
            for (int c=0; c<n; c++)
                std::cout << C[c] << "\t";
            std::cout << obj << std::endl;
            for (int r=0; r<m; r++) {
                for (int c=0; c<n; c++)
                    std::cout << A[r][c] << "\t";
                std::cout << B[r] << std::endl;
            }
            std::cout << std::endl;
        }
    }
    return false;
}

/*##########################################################
 The main functions for solving using the simplex algorithm
 simplex Phase I : for finding an initial basis
 simplex Phase II: for solving the problem
 ##########################################################*/
void LP_solver::set_initial_basis(){
    
    mat_d & A = LP.A;
    //vec_d & b = LP.b;// is assumed having >= entries
    size_t num_rows = A.size();
    size_t num_cols = A[0].size();
    size_t basis_size = basis.size();
    for (size_t j=0; j<num_cols; j++) { //First sweep: Find the initial basic variables.
        int r = is_unit_column (A, j);
        if (r >= 0 && basis[r] < 0) {
            basis[r] = (int)j;
            basis_size++;
        }
    }
    // if basis_size ==num_rows, no need to go further,
    if (basis_size == num_rows ) { return;}
    // otherwise we add artificial variables by adding columns to matrix A
    // this operation well be performed in the phase I of simplex method
}

// Returns: -1 if unbounded; 0 if infeasible; 1 if feasible & bounded.
void LP_solver::run_simplex_phase_I(){
    size_t m = LP.A.size();
    size_t n = LP.A[0].size();
    
    if (m == basis.size()){
        //  is_unbounded = false;
        is_unfeasible = false;
        return;
    }
    
    // we use notation A1, b1, c1 : for Phase I.
    mat_d& A1 =LP.A;
    vec_d& b1 =LP.b;
    vec_d c1 (n, 0);    // new objective vector for phase I
    
    
    size_t num_artificial = m - basis.size();
    
    if (num_artificial > 0){
        int j = 0;
        size_t i = 0;
        for (vec_d& curr_row : A1){
            if (basis[i] < 0){
                concatenate_kronecker(curr_row, num_artificial, j);
                basis[i] =(int) n+j;
                j++;
            }
            else
                concatenate_zeros(curr_row, num_artificial );
            i++;
        }
    }
    
    size_t n1 = n + num_artificial ;        // Phase I need artificial variables !!
    c1.resize (n1, 1);   // Adjust size of objective vector adding artificial variables
    optimum = 0.;
    
    for (int i=0; i<m; i++) {LP.b[i]+=EPS;}
    canonicalize (A1, b1, c1, basis.basis, optimum);  // convert to canonical form
    is_unbounded = pivoting (A1, b1, c1, basis.basis, optimum, 1000);  // pivoting
    is_unfeasible = (fabs(optimum) < EPS) ?  false: true;
    // Initial b.f.s. found, prepare for Phase II.
    for (vec_d& curr_row : A1)
        curr_row.resize(n);
}

void LP_solver::run_simplex_phase_II(){
    size_t m = LP.A.size();
    // basis is already set thanks to Phase I of simplex !!
    // we use notation A2, b2, c2 : for Phase II.
    mat_d& A2 =LP.A;
    vec_d& b2 =LP.b;
    vec_d& c2 =LP.c;
    // Phase II.
    canonicalize (A2, b2, c2, basis.basis, optimum);
    is_unbounded = pivoting (A2, b2, c2, basis.basis, optimum, 1000);
    if(!is_unbounded) is_solved =true;
    std::for_each(solution.begin(),solution.end(),[](double &x){x=0.0;});
    for (int r=0; r<m; r++)              // r.h.s. is the basic solution
        solution[basis[r]] = b2[r];
}
}
#endif /* SLE_simplex_h */
