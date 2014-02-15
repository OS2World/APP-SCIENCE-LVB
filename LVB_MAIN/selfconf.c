/* LVB (c) Copyright 2003-2006 by Daniel Barker.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** selfconf.c - self-configuration functions ********** */

/* Define arrays (and their extents) used in the search for best simulated
 * annealing parameters */

#include "lvb.h"

static const char *rcsid
 = " $Id: selfconf.c,v 1.28 2006/02/06 19:55:46 db60 Exp $";

#define MAXACCEPT_5_LIM ((5 + 50) / 2)

#define T0_CNT 21
/* members of the array must be in decreasing order */
static const double t0_arr[T0_CNT] =
{
    1,			/* 10^ 0 */
    0.630957,		/* 10^-0.2 */
    0.398107,		/* 10^-0.4 */
    0.251189,		/* 10^-0.6 */
    0.158489,		/* 10^-0.8 */
    0.1,		/* 10^-1 */
    0.0630957,		/* 10^-1.2 */
    0.0398107,		/* 10^-1.4 */
    0.0251189,		/* 10^-1.6 */
    0.0158489,		/* 10^-1.8 */
    0.01,		/* 10^-2 */
    0.00630957,		/* 10^-2.2 */
    0.00398107,		/* 10^-2.4 */
    0.00251189,		/* 10^-2.6 */
    0.00158489,		/* 10^-2.8 */
    0.001,		/* 10^-3 */
    0.000630957,	/* 10^-3.2 */
    0.000398107,	/* 10^-3.4 */
    0.000251189,	/* 10^-3.6 */
    0.000158489,	/* 10^-3.8 */
    0.0001		/* 10^-4 */
};

#define T1_TO_T0_CNT 16
/* members of the array must be in decreasing order */
static const double t1_to_t0_arr[T1_TO_T0_CNT] =
{
    0.99,
    0.500,
    0.099,
    0.0500,
    0.0400,
    0.0300,
    0.0200,
    0.0160,
    0.0150,
    0.0140,
    0.0130,
    0.0120,
    0.0110,
    0.0100,
    0.00995,
    0.0099
};

#define MAXPROPOSE_CNT 15
/* members of the array must be in increasing order */
static const long maxpropose_arr[MAXPROPOSE_CNT] =
{
    50,
    100,
    200,
    300,
    400,
    500,
    600,
    700,
    800,
    900,
    1000,
    2000,
    3000,
    4000,
    5000
};

#define MAXACCEPT_CNT 2
/* members of the array must be in increasing order */
static const long maxaccept_arr[MAXACCEPT_CNT] = { 5, 50 };

#define MAXFAIL_CNT 15
/* members of the array must be in increasing order */
static const long maxfail_arr[MAXFAIL_CNT] =
{
     1,
     2,
     3,
     4,
     5,
     6,
     7,
     8,
     9,
    10,
    20,
    40,
    60,
    80,
    100
};

#define RUNS_CNT 10
/* members of the array must be in increasing order */
static const long runs_arr[RUNS_CNT] =
{
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10
};

typedef enum
{
    STEPS,
    TREES
} Response_type;

static double estimate_response(long maxaccept, double a, double b,
 double c, double d)
/*@modifies nothing@*/
/* Calculate response variable in a fitted equation for factor
 * maxaccept, given a, b, c and d. For explanationsof a, b, c and d see
 * Table 4.4 of D. Barker (1999), PhD Thesis, University of
 * Edinburgh */
{
    double resp;	/* fitted value */

    /* To obtain fitted equations, in the thesis, only two levels were
     * used for maxaccept. We decide which is most appropriate here,
     * and then, for the current function, pretend maxaccept always had
     * that value. */
    if (maxaccept <= MAXACCEPT_5_LIM)
        maxaccept = 5;
    else
    	maxaccept = 50;

    if (maxaccept == 5)
    {
        resp = a + b + c + d;
    }
    else
    {
        resp = a + b - c - d;
    }

    return resp;

}	/* end estimate_response() */

static void abcd_lnsteps(double ln_t0, double ln_t1_to_t0,
 double ln_maxpropose, double ln_maxfail, long runs,
 /*@out@*/ double *a_lnsteps_p, /*@out@*/ double *b_lnsteps_p,
 /*@out@*/ double *c_lnsteps_p, /*@out@*/ double *d_lnsteps_p)
/*@modifies *a_lnsteps_p, *b_lnsteps_p, *c_lnsteps_p, *d_lnsteps_p@*/
/* Calculate a_lnSteps, b_lnSteps, c_lnSteps and d_lnSteps for the given
 * values of ln(t0), ln(t1/t0), ln(maxpropose), ln(maxfail) and runs, passed
 * as ln_t0, ln_t1_to_t0, ln_maxpropose, ln_maxfail and runs respectively.
 * On return, a_lnSteps, b_lnSteps, c_lnSteps and d_lnSteps are left in
 * *a_lnsteps_p, *b_lnsteps_p, *c_lnsteps_p and *d_lnsteps_p respectively. */
{
    double a_lnsteps;	/* a_lnSteps */
    double b_lnsteps;	/* b_lnSteps */
    double c_lnsteps;	/* c_lnSteps */
    double d_lnsteps;	/* d_lnSteps */

    /* calculate the required numbers */

    a_lnsteps = 8.77390
     - 0.003750 * runs
     - 0.01829 * ln_t1_to_t0
     - 0.003342 * ln_maxpropose
     - 0.01106 * ln_maxfail
     - 0.017079 * ln_t0
     - 0.002349 * ln_t0 * ln_t1_to_t0
     + 0.006575 * ln_t0 * ln_maxpropose
     + 0.007075 * ln_t0 * ln_maxfail
     + 0.009346 * ln_t1_to_t0 * ln_maxpropose
     + 0.014286 * ln_t1_to_t0 * ln_maxfail
     - 0.001108 * ln_maxpropose * ln_maxfail
     + 0.000863 * ln_t0 * ln_t1_to_t0 * ln_maxpropose
     + 0.001250 * ln_t0 * ln_t1_to_t0 * ln_maxfail
     - 0.000709 * ln_t0 * ln_maxpropose * ln_maxfail
     - 0.000343 * ln_t1_to_t0 * ln_maxpropose * ln_maxfail
     + 0.000039 * ln_t0 * ln_t1_to_t0
     * ln_maxpropose * ln_maxfail;

    b_lnsteps = 0.95173 - 0.03933 * ln_t1_to_t0
     + 0.014908 * ln_maxpropose
     + 0.03144 * ln_maxfail
     - 0.017762 * ln_t0
     - 0.005040 * ln_t0 * ln_t1_to_t0
     + 0.003745 * ln_t1_to_t0 * ln_maxfail
     + 0.002415 * ln_t0 * ln_maxpropose
     + 0.001022 * ln_t0 * ln_maxfail
     - 0.004006 * ln_maxpropose * ln_maxfail
     + 0.000134 * ln_t0 * ln_maxpropose * ln_maxfail
     + 0.005693 * ln_t1_to_t0 * ln_maxpropose
     + 0.000658 * ln_t0 * ln_t1_to_t0 * ln_maxpropose
     + 0.000550 * ln_t0 * ln_t1_to_t0 * ln_maxfail
     + 0.000025 * ln_t1_to_t0
     * ln_maxpropose * ln_maxfail;
    
    c_lnsteps = 0.05133 + 0.00415 * ln_t1_to_t0
     - 0.009604 * ln_maxpropose
     - 0.01421 * ln_maxfail
     + 0.004311 * ln_t1_to_t0 * ln_maxpropose
     + 0.008893 * ln_t0
     - 0.002896 * ln_t0 * ln_t1_to_t0
     + 0.006561 * ln_t1_to_t0 * ln_maxfail
     - 0.000586 * ln_t0 * ln_maxpropose
     + 0.001012 * ln_t0 * ln_t1_to_t0 * ln_maxpropose
     - 0.000760 * ln_t0 * ln_maxfail
     + 0.001254 * ln_t0 * ln_t1_to_t0 * ln_maxfail
     + 0.001863 * ln_maxpropose * ln_maxfail
     - 0.002669 * ln_t1_to_t0 * ln_maxpropose * ln_maxfail
     + 0.000069 * ln_t0 * ln_maxpropose * ln_maxfail
     - 0.000297 * ln_t0 * ln_t1_to_t0
     * ln_maxpropose * ln_maxfail;
    
    d_lnsteps = 0.03721 - 0.002316 * ln_t0
     + 0.00283 * ln_t1_to_t0
     - 0.001547 * ln_t0 * ln_t1_to_t0
     - 0.010373 * ln_maxfail
     + 0.003334 * ln_t1_to_t0 * ln_maxfail
     - 0.006370 * ln_maxpropose
     + 0.000658 * ln_t0 * ln_maxpropose
     - 0.000753 * ln_t1_to_t0 * ln_maxpropose
     + 0.000286 * ln_t0 * ln_t1_to_t0 * ln_maxpropose
     + 0.000694 * ln_t0 * ln_maxfail
     + 0.000244 * ln_t0 * ln_t1_to_t0 * ln_maxfail
     + 0.001743 * ln_maxpropose * ln_maxfail
     - 0.000673 * ln_t1_to_t0
     * ln_maxpropose * ln_maxfail;

    /* make results available to caller */
    *a_lnsteps_p = a_lnsteps;
    *b_lnsteps_p = b_lnsteps;
    *c_lnsteps_p = c_lnsteps;
    *d_lnsteps_p = d_lnsteps;

}	/* end abcd_lnsteps() */

static void abcd_lntrees(double ln_t0, double ln_t1_to_t0,
 double ln_maxpropose, double ln_maxfail, long runs,
 /*@out@*/ double *a_lntrees_p, /*@out@*/ double *b_lntrees_p,
 /*@out@*/ double *c_lntrees_p, /*@out@*/ double *d_lntrees_p)
/*@modifies *a_lntrees_p, *b_lntrees_p, *c_lntrees_p, *d_lntrees_p@*/
/* Calculate a_lnTrees, b_lnTrees, c_lnTrees and d_lnTrees for the given
 * values of ln(t0), ln(t1/t0), ln(maxpropose), ln(maxfail) and runs, passed
 * as ln_t0, ln_t1_to_t0, ln_maxpropose and ln_maxfail and runs respectively.
 * On return, a_lnTrees, b_lnTrees, c_lnTrees and d_lnTrees are left in
 * *a_lnsteps_p, *b_lnsteps_p, *c_lnsteps_p and *d_lnsteps_p respectively. */
{
    double a_lntrees;	/* a_lnTrees */
    double b_lntrees;	/* b_lnTrees */
    double c_lntrees;	/* c_lnTrees */
    double d_lntrees;	/* d_lnTrees */

    /* calculate the required numbers */

    a_lntrees = -0.3007
     + 0.256130 * runs
     - 0.07438 * ln_t0
     - 0.24180 * ln_t1_to_t0
     + 1.03934 * ln_maxpropose
     + 1.09576 * ln_maxfail
     - 0.001081 * ln_t0 * ln_t1_to_t0
     + 0.012133 * ln_t0 * ln_maxfail
     + 0.011791 * ln_t1_to_t0 * ln_maxpropose
     - 0.00180 * ln_t1_to_t0 * ln_maxfail
     - 0.008501 * ln_maxpropose * ln_maxfail
     + 0.009295 * ln_t1_to_t0 * ln_maxpropose * ln_maxfail
     - 0.006515 * ln_t0 * ln_maxpropose
     - 0.002036 * ln_t0 * ln_t1_to_t0 * ln_maxpropose
     - 0.001389 * ln_t0 * ln_t1_to_t0 * ln_maxfail
     + 0.001517 * ln_t0 * ln_maxpropose * ln_maxfail
     + 0.000396 * ln_t0 * ln_t1_to_t0
     * ln_maxpropose * ln_maxfail;

    b_lntrees = -0.0809
     + 0.03292 * ln_t0
     - 0.00803 * ln_t1_to_t0
     + 0.00477 * ln_maxpropose
     + 0.01127 * ln_maxfail
     + 0.019227 * ln_t0 * ln_t1_to_t0
     - 0.011291 * ln_t0 * ln_maxpropose
     - 0.05312 * ln_t1_to_t0 * ln_maxfail
     - 0.001976 * ln_t1_to_t0 * ln_maxpropose
     - 0.003102 * ln_t0 * ln_t1_to_t0 * ln_maxpropose
     + 0.001011 * ln_maxpropose * ln_maxfail
     + 0.011372 * ln_t1_to_t0 * ln_maxpropose * ln_maxfail
     - 0.011838 * ln_t0 * ln_maxfail
     - 0.005459 * ln_t0 * ln_t1_to_t0 * ln_maxfail
     + 0.002257 * ln_t0 * ln_maxpropose * ln_maxfail
     + 0.000544 * ln_t0 * ln_t1_to_t0
     * ln_maxpropose * ln_maxfail;
    
    c_lntrees = -0.4297
     - 0.13979 * ln_t0
     - 0.28009 * ln_t1_to_t0
     + 0.08720 * ln_maxpropose
     + 0.14833 * ln_maxfail
     - 0.021647 * ln_maxpropose * ln_maxfail
     + 0.012504 * ln_t0 * ln_maxpropose
     + 0.026390 * ln_t0 * ln_maxfail
     + 0.021015 * ln_t1_to_t0 * ln_maxpropose
     + 0.001070 * ln_t0 * ln_t1_to_t0
     - 0.00114 * ln_t1_to_t0 * ln_maxfail
     + 0.009566 * ln_t1_to_t0 * ln_maxpropose * ln_maxfail
     - 0.002718 * ln_t0 * ln_t1_to_t0 * ln_maxpropose
     - 0.003037 * ln_t0 * ln_t1_to_t0 * ln_maxfail
     - 0.002744 * ln_t0 * ln_maxpropose * ln_maxfail
     + 0.000855 * ln_t0 * ln_t1_to_t0
     * ln_maxpropose * ln_maxfail;
    
    d_lntrees = -0.2173
     - 0.00280 * ln_t0
     - 0.03429 * ln_t1_to_t0
     + 0.022414 * ln_t0 * ln_t1_to_t0
     + 0.04485 * ln_maxfail
     - 0.05170 * ln_t1_to_t0 * ln_maxfail
     + 0.03971 * ln_maxpropose
     - 0.001530 * ln_t0 * ln_maxpropose
     + 0.004010 * ln_t1_to_t0 * ln_maxpropose
     - 0.004063 * ln_t0 * ln_t1_to_t0 * ln_maxpropose
     - 0.005449 * ln_t0 * ln_maxfail
     - 0.006983 * ln_t0 * ln_t1_to_t0 * ln_maxfail
     - 0.007646 * ln_maxpropose * ln_maxfail
     + 0.011459 * ln_t1_to_t0 * ln_maxpropose * ln_maxfail
     + 0.000412 * ln_t0 * ln_maxpropose * ln_maxfail
     + 0.000976 * ln_t0 * ln_t1_to_t0
     * ln_maxpropose * ln_maxfail;
    
    /* make results available to caller */
    *a_lntrees_p = a_lntrees;
    *b_lntrees_p = b_lntrees;
    *c_lntrees_p = c_lntrees;
    *d_lntrees_p = d_lntrees;
 
}	/* end abcd_lntrees() */

static double get_fitted_value(double t0, double t1, long maxpropose,
 long maxaccept, long maxfail, long runs, Response_type resp)
/*@*/
{
    double a_lnresp;		/* a_lnSteps or a_lnTrees */
    double b_lnresp;		/* b_lnSteps or b_lnTrees */
    double c_lnresp;		/* c_lnSteps or c_lnTrees */
    double d_lnresp;		/* d_lnSteps or d_lnTrees */
    double lnresp;		/* estimate for ln Steps or ln Trees */
    double ln_t0; 		/* ln T0 */
    double ln_t1_to_t0;		/* ln T1 */
    double ln_maxpropose;	/* ln maxpropose */
    double ln_maxfail;		/* ln maxfail */
    double val;			/* return value */

    ln_t0 = log_wrapper(t0);
    ln_t1_to_t0 = log_wrapper(t0 / t1);
    ln_maxpropose = log_wrapper((double) maxpropose);
    ln_maxfail = log_wrapper((double) maxfail);

    lvb_assert((resp == TREES) || (resp == STEPS));
    if (resp == STEPS)
    {
        abcd_lnsteps(ln_t0, ln_t1_to_t0, ln_maxpropose, ln_maxfail,
         runs, &a_lnresp, &b_lnresp, &c_lnresp, &d_lnresp);
        lnresp = estimate_response(maxaccept, a_lnresp, b_lnresp,
	 c_lnresp, d_lnresp);
        val = exp_wrapper(lnresp);
    }
    else
    {
        abcd_lntrees(ln_t0, ln_t1_to_t0, ln_maxpropose, ln_maxfail,
         runs, &a_lnresp, &b_lnresp, &c_lnresp, &d_lnresp);
        lnresp = estimate_response(maxaccept, a_lnresp, b_lnresp,
	 c_lnresp, d_lnresp);
        val = exp_wrapper(lnresp);
    }
    return val;

}	/* end get_fitted_value() */

double get_predicted_length(double t0, double t1,
 long maxpropose, long maxaccept, long maxfail, long runs)
/*@*/
/* Return the predicted best tree length for a large matrix with the
 * simulated annealing parameters given by t0, t1, maxpropose,
 * maxaccept, maxfail and runs. */
{
   double val;	/* return value */

   val = get_fitted_value(t0, t1, maxpropose, maxaccept, maxfail,
    runs, STEPS);
   return val;

}	/* end get_predicted_length() */

double get_predicted_trees(double t0, double t1,
 long maxpropose, long maxaccept, long maxfail, long runs)
/*@*/
/* Return the predicted number of trees to be examined
 * with the simulated annealing parameters given by t0, t1, maxpropose,
 * maxaccept, maxfail and runs. */
{
   double val;	/* return value */

   val = get_fitted_value(t0, t1, maxpropose, maxaccept, maxfail,
    runs, TREES);
   return val;

}	/* end get_predicted_trees() */

Anneal_params_struct self_configure(long max_trees)
/*@globals t0_arr, t1_to_t0_arr, maxpropose_arr, maxaccept_arr, maxfail_arr,
 runs_arr@*/ /*@modifies nothing@*/
/* For an analysis predicted to involve
 * examination of at most max_trees trees, return simulated annealing
 * parameters predicted to
 * find the shortest tree. Also return the associated details
 * mat_type, predicted best tree length and predicted number of trees
 * examined (which may be less than max_trees). */
{
    Anneal_params_struct params;	/* current set of parameters */
    Anneal_params_struct best_params;	/* best set of parameters */
    long i;				/* loop counter */
    long j;				/* loop counter */
    long k;				/* loop counter */
    long m;				/* loop counter */
    long n;				/* loop counter */
    long p;				/* loop counter */

    best_params.predicted_length = LONG_MAX;

    /* Set annealing parameters in best_params to the values that will
     * give the quickest possible run. Without this, if the user asks
     * for an astonishingly quick analysis, best_params will not get
     * set in the loops below (due to the impossibility of the
     * request). */
    best_params.t0 = t0_arr[T0_CNT - 1];
    best_params.t1 = best_params.t0 * t1_to_t0_arr[T1_TO_T0_CNT - 1];
    best_params.maxpropose = maxpropose_arr[0];
    best_params.maxaccept = maxaccept_arr[0];
    best_params.maxfail = maxfail_arr[0];
    best_params.runs = runs_arr[0];

    for (i = 0; i < T0_CNT; i++)
    {
        params.t0 = t0_arr[i];
        for (j = 0; j < T1_TO_T0_CNT; j++)
	{
	    params.t1 = params.t0 * t1_to_t0_arr[j];
	    for (k = 0; k < MAXPROPOSE_CNT; k++)
	    {
	        params.maxpropose = maxpropose_arr[k];
	        for (m = 0; m < MAXACCEPT_CNT; m++)
		{
		    params.maxaccept = maxaccept_arr[m];
	            for (n = 0; n < MAXFAIL_CNT; n++)
                    {
		        params.maxfail = maxfail_arr[n];
		        for (p = 0; p < RUNS_CNT; p++)
		        {
			    params.runs = runs_arr[p];
			    params.predicted_length =
			     get_predicted_length(params.t0, params.t1,
                             params.maxpropose,
			     params.maxaccept, params.maxfail,
			     params.runs);
			    if (params.predicted_length
			     < best_params.predicted_length)
			    {
			        params.predicted_trees =
				 get_predicted_trees(params.t0, params.t1,
                                 params.maxpropose,
			         params.maxaccept, params.maxfail,
			         params.runs);
				if (params.predicted_trees < (double) max_trees)
				{
				    best_params = params;
				}
			    }
			}
		    }
		}
	    }
	}
    }

    /* If best_params has been unaffected by the above loops, we still
     * have to set its fields that show predicted tree length and run
     * time. */
    if (best_params.predicted_length == LONG_MAX)
    {
	best_params.predicted_length = get_predicted_length(best_params.t0,
         best_params.t1, best_params.maxpropose,
	 best_params.maxaccept, best_params.maxfail, best_params.runs);
	best_params.predicted_trees = get_predicted_trees(best_params.t0,
         best_params.t1, best_params.maxpropose,
	 best_params.maxaccept, best_params.maxfail, best_params.runs);
    }

    return best_params;

}	/* end self_configure() */
