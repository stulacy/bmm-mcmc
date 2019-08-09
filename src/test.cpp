#include <Rcpp.h>
#include <stdio.h>
#include "ls_source/lp_lib.h"

using namespace Rcpp;

void lp_transbig_edit(int, int, double *, double *);

// [[Rcpp::export]]
IntegerMatrix lpsolve(NumericMatrix x) {
    int rcount = x.nrow();
    int ccount = x.ncol();
    double objective[1+(rcount * ccount)];
    double solution[rcount * ccount];
    
    objective[0] = 0;
    for (int i = 0; i < rcount; ++i) {
        for (int j = 0; j < rcount; ++j) {
            objective[(i * ccount) + j + 1] = x(i, j);
            solution[(i * ccount) + j] = 0;
        }
    }
    
    lp_transbig_edit(
                rcount, 
                ccount, 
                objective, 
                solution);
    
    // Convert into matrix
    IntegerMatrix sol(rcount, ccount);
    for (int i = 0; i < rcount; ++i) {
        for (int j = 0; j < ccount; ++j) {
            sol(i, j) = solution[(i*ccount) + j];
        }
    }
        
    return sol;
}   

void lp_transbig_edit (
              int r_count,           /* Number of rows             */
              int c_count,           /* Number of columns          */
              double *costs,                  /* Objective function         */
              double *solution)               /* Result of call             */
{
    long i;              /* Iteration variable       */
    long result;         /* Holds result of calls    */
    long this_element;   /* Which are we looking at? */
    lprec *lp;           /* Structure to hold the lp */
    double *row_vals;    /* Holds the values for row-type constraints */
    int *col_inds;       /* Holds locations for col-type constraints  */
    double *col_vals;    /* Holds the values for col-type constraints */
    int *row_inds;       /* Holds locations for row-type constraints  */
    
    long col_ind_ctr, row_ind_ctr;
    long rc = r_count, cc = c_count;
    
    /*
    ** Make an empty lp with r_count x c_count variables. If it fails, return.
    */
    lp = make_lp ((int) 0, r_count * c_count);
    
    if (lp == (lprec *) NULL)
        return;
    
    set_verbose (lp, 1); /* CRITICAL */
    
    /*
    ** "Costs" is already a vector. Set the objective function. Return on fail.
    */
    result = set_obj_fn (lp, costs);
    if (result == 0)
        return;
    
    // Minimising loss
    set_minim (lp);
    
    // Add constraints
    row_vals = (double *) calloc (cc, sizeof (double));
    col_inds = (int *) calloc (cc, sizeof (int));
    
    int rsign = 3;
    int csign = 3;
    int rconstraint = 1;
    int cconstraint = 1;
    
    for (row_ind_ctr = 0L; row_ind_ctr < rc; row_ind_ctr++)
    {
        for (col_ind_ctr = 0; col_ind_ctr < cc; col_ind_ctr++) {
            row_vals[col_ind_ctr] = 1.0;
            this_element = 1 + (col_ind_ctr * rc) + row_ind_ctr;
            col_inds[col_ind_ctr] = this_element;
        }
        add_constraintex (lp, cc, row_vals, col_inds, rsign, rconstraint);
    }
    
    free (row_vals);
    free (col_inds);
    
    col_vals = (double *) calloc (rc, sizeof (double));
    row_inds = (int *) calloc (rc, sizeof (int));
    
    for (col_ind_ctr = 0L; col_ind_ctr < cc; col_ind_ctr++)
    {
        for (row_ind_ctr = 0; row_ind_ctr < rc; row_ind_ctr++) {
            col_vals[row_ind_ctr] = 1.0;
            this_element = 1 + row_ind_ctr + col_ind_ctr * rc;
            row_inds[row_ind_ctr] = this_element;
        }
        add_constraintex (lp, rc, col_vals, row_inds, csign, cconstraint);
    }
    free (col_vals);
    free (row_inds);
    
    set_add_rowmode (lp, FALSE);
    
    /*
    ** Set integers.
    */
    for (i = 1; i <= (r_count * c_count); ++i)
        set_int (lp, i, 1); /* Variable in ith element of integers */
    
    if (solve(lp) != 0) {
        return;
    }
    get_variables(lp, solution);
    delete_lp (lp);
}

