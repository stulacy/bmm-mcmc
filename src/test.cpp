//// [[Rcpp::depends(RcppArmadillo)]]
//
//#include <RcppArmadilloExtensions/sample.h>
//#include <stdio.h>
//#include "lp_lib.h"
//
//using namespace Rcpp;
//
//// This is a simple example of exporting a C++ function to R. You can
//// source this function into an R session using the Rcpp::sourceCpp
//// function (or via the Source button on the editor toolbar). Learn
//// more about Rcpp at:
////
////   http://www.rcpp.org/
////   http://adv-r.had.co.nz/Rcpp.html
////   http://gallery.rcpp.org/
////
//void lp_transbig_edit(int, int, double *, double *, int *, double *);
//
//// [[Rcpp::export]]
//NumericVector timesTwo(NumericVector x) {
//    arma::mat foo(2,2);
//    foo(0, 0) = -1.267420;
//    foo(0, 1) = -1.267420;
//    foo(1, 0) = 2.074419;
//    foo(1, 1) = 2.074419;
//    Rcout << "Foo: " << foo  <<         "\n";
//    
//    int rcount = 2;
//    int ccount = 2;
//    double objective[5] = {0, -1.267420, -1.267420, 2.074419, 2.074419};
//    double *objval = 0;
//    int integers[4] = {0, 0, 0, 0};
//    double solution[4] = {0.0, 0.0, 0.0, 0.0};
//    
//    lp_transbig_edit(
//                rcount, 
//                ccount, 
//                objective, 
//                objval,
//                integers, 
//                solution);
//    return x * 2;
//}   
//
//void lp_transbig_edit (
//              int r_count,           /* Number of rows             */
//              int c_count,           /* Number of columns          */
//              double *costs,                  /* Objective function         */
//              double *obj_val,                /* Objective function value   */
//              int *integers,          /* Which vars. are integer?   */
//              double *solution)               /* Result of call             */
//{
//    long i;              /* Iteration variable       */
//    long result;         /* Holds result of calls    */
//    long this_element;   /* Which are we looking at? */
//    lprec *lp;           /* Structure to hold the lp */
//    double *row_vals;    /* Holds the values for row-type constraints */
//    int *col_inds;       /* Holds locations for col-type constraints  */
//    double *col_vals;    /* Holds the values for col-type constraints */
//    int *row_inds;       /* Holds locations for row-type constraints  */
//    
//    long col_ind_ctr, row_ind_ctr;
//    long rc = r_count, cc = c_count;
//    
//    /*
//    ** Make an empty lp with r_count x c_count variables. If it fails, return.
//    */
//    lp = make_lp ((int) 0, r_count * c_count);
//    
//    if (lp == (lprec *) NULL)
//        return;
//    
//    set_verbose (lp, 1); /* CRITICAL */
//    
//    /*
//    ** "Costs" is already a vector. Set the objective function. Return on fail.
//    */
//    result = set_obj_fn (lp, costs);
//    if (result == 0)
//        return;
//    
//    // Minimising loss
//    set_minim (lp);
//    
//    /*
//    ** Add constraints. There are r_count row-type constraints, plus c_count
//    ** col_type constraints.
//    */
//    row_vals = (double *) calloc (cc, sizeof (double));
//    col_inds = (int *) calloc (cc, sizeof (int));
//    
//    int rsign = 3;
//    int csign = 3;
//    int rconstraint = 1;
//    int cconstraint = 1;
//    
//    for (row_ind_ctr = 0L; row_ind_ctr < rc; row_ind_ctr++)
//    {
//        for (col_ind_ctr = 0; col_ind_ctr < cc; col_ind_ctr++) {
//            row_vals[col_ind_ctr] = 1.0;
//            this_element = 1 + (col_ind_ctr * rc) + row_ind_ctr;
//            col_inds[col_ind_ctr] = this_element;
//        }
//        add_constraintex (lp, cc, row_vals, col_inds, rsign, rconstraint);
//    }
//    
//    free (row_vals);
//    free (col_inds);
//    
//    col_vals = (double *) calloc (rc, sizeof (double));
//    row_inds = (int *) calloc (rc, sizeof (int));
//    
//    for (col_ind_ctr = 0L; col_ind_ctr < cc; col_ind_ctr++)
//    {
//        for (row_ind_ctr = 0; row_ind_ctr < rc; row_ind_ctr++) {
//            col_vals[row_ind_ctr] = 1.0;
//            this_element = 1 + row_ind_ctr + col_ind_ctr * rc;
//            row_inds[row_ind_ctr] = this_element;
//        }
//        add_constraintex (lp, rc, col_vals, row_inds, csign, cconstraint);
//    }
//    free (col_vals);
//    free (row_inds);
//    
//    set_add_rowmode (lp, FALSE);
//    
//    /*
//    ** Set integers.
//    */
//    for (i = 0; i < r_count * c_count; i++)
//        set_int (lp, integers[i], 1); /* Variable in ith element of integers */
//    
//    if (solve(lp) != 0) {
//        return;
//    }
//    
//    *obj_val = get_objective (lp);
//    get_variables (lp, solution);
//    delete_lp (lp);
//}
//
//