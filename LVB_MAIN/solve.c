/* LVB (c) Copyright 2003-2006 by Daniel Barker.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** solve.c - solving functions ********** */

#include "lvb.h"

static const char *rcsid = "$Id: solve.c,v 1.47 2006/02/06 19:55:47 db60 Exp $";

static void lenlog(FILE *lengthfp, long iteration, long length)
/* write a message to file pointer lengthfp; iteration gives current iteration;
 * crash verbosely on write error */
{
    fprintf(lengthfp, "%-15ld%ld\n", iteration, length); 
    if (ferror(lengthfp))
    {
	crash("file error when logging search progress");
    }

} /* end lenlog() */

long anneal(Treestack *bstackp, const Branch *const inittree, long root,
 const double t0, const double t1, const long maxaccept, const long maxpropose,
 const long maxfail, FILE *const lenfp, unsigned char **bmat, long m, long n,
 const long *weights, long *current_iter, Lvb_bool log_progress)
/* seek parsimonious tree from initial tree in inittree (of root root)
 * with initial temperature t0, and subsequent temperatures obtained by
 * multiplying the current temperature by (t1 / t0) ** n * t0 where n is
 * the ordinal number of this temperature, after at least maxaccept changes
 * have been accepted or maxpropose changes have been proposed, whichever is
 * sooner;
 * return the length of the best tree(s) found after maxfail consecutive
 * temperatures have led to no new accepted solution;
 * lenfp is for output of current tree length and associated details;
 * *current_iter should give the iteration number at the start of this call and
 * will be used in any statistics sent to lenfp, and will be updated on
 * return */
{
    long accepted = 0;		/* changes accespted */
    Lvb_bool dect;		/* should decrease temperature */
    double deltah;		/* change in energy (1 - C.I.) */
    long deltalen;		/* change in length with new tree */
    long failedcnt = 0; 	/* "failed count" for temperatures */
    long iter = 0;		/* iteration of mutate/evaluate loop */
    long len;			/* length of current tree */
    long prev_len = UNSET;	/* length of previous tree */
    long lenbest;		/* bet length found so far */
    long lendash;		/* length of proposed new tree */
    long lenmin;		/* minimum length for any tree */
    double ln_t;		/* ln(current temperature) */
    long t_n = 0;		/* ordinal number of current temperature */
    Lvb_bool newtree;		/* accepted a new configuration */
    double pacc;		/* prob. of accepting new config. */
    Lvb_bool probaccd;		/* have accepted based on Pacc */
    long proposed = 0;		/* trees proposed */
    double r_lenmin;		/* minimum length for any tree */
    long rootdash;		/* root of new configuration */
    double t = t0;		/* current temperature */
    double t1_to_t0;		/* T1:T0 ratio */
    Branch *x;			/* current configuration */
    Branch *xdash;		/* proposed new configuration */
    extern Dataptr matrix;	/* data matrix */

    /* "local" dynamic heap memory */
    x = treealloc(brcnt(matrix->n), matrix->m);
    xdash = treealloc(brcnt(matrix->n), matrix->m);

    treecopy(x, inittree);	/* current configuration */
    len = getplen(x, root, m, n, weights);
    dect = LVB_FALSE;		/* made LVB_TRUE as necessary at end of loop */

    /* if t1_to_t0 is going to be very small, set it to zero without
     * calculating it, to avoid underflow. Use this relation:
     *     if T1 / T0 < eps,
     *     ln (T1/T0) < ln eps
     * i.e.,
     *     ln T1 - ln T0 < ln eps */
    lvb_assert ((t0 >= LVB_EPS) && (t1 >= LVB_EPS) && (t1 <= t0)
     && (t0 <= 1.0));
    if (log_wrapper(t1) - log_wrapper(t0) < log_wrapper(LVB_EPS))
	t1_to_t0 = 0.0;
    else
	t1_to_t0 = t1 / t0;

    lenbest = len;
    treestack_push(bstackp, inittree, root);	/* init. tree initially best */
    if ((log_progress == LVB_TRUE) && (*current_iter == 0))
    {
        fprintf(lenfp, "\nRearrangement: Length:\n");
    }

    lenmin = getminlen(matrix);
    r_lenmin = (double) lenmin;

    while (1)
    {
        if ((log_progress == LVB_TRUE) && ((len != prev_len)
	    || ((*current_iter % STAT_LOG_INTERVAL) == 0)))
        {
	    lenlog(lenfp, *current_iter, len);
        }
	prev_len = len;
	*current_iter += 1;

	/* occasionally re-root, to prevent influence from root position */
	if ((*current_iter % REROOT_INTERVAL) == 0)
	    root = arbreroot(x, root);

	lvb_assert(t > DBL_EPSILON);
	newtree = LVB_FALSE;
	probaccd = LVB_FALSE;

	/* mutation: alternate between the two mutation functions */
	if (iter % 2)
	{
	    rootdash = root;
	    mutate_spr(xdash, x, root);	/* global change */
	}
	else
	{
	    rootdash = root;
	    mutate_nni(xdash, x, root);	/* local change */
	}

	lendash = getplen(xdash, rootdash, m, n, weights);
	lvb_assert (lendash >= 1L);
	deltalen = lendash - len;
	deltah = (r_lenmin / (double) len) - (r_lenmin / (double) lendash);
	if (deltah > 1.0)	/* getminlen() problem with ambiguous sites */
	    deltah = 1.0;
	if (deltalen <= 0)	/* accept the change */
	{
	    if (lendash <= lenbest)	/* store tree if new */
	    {
		if (lendash < lenbest)	/* very best so far */
		treestack_clear(bstackp);	/* discard old bests */
		if (treestack_push(bstackp, xdash, rootdash) == 1)
		newtree = LVB_TRUE;	/* new */
	    }
	    /* update current tree and its stats */
	    prev_len = len;
	    len = lendash;
	    treeswap(&x, &root, &xdash, &rootdash);

	    if (lendash < lenbest)	/* very best so far */
		lenbest = lendash;
	}
	else	/* poss. accept change for the worse */
	{
	    /* Mathematically,
	     *     Pacc = e ** (-1/T * deltaH)
	     *     therefore ln Pacc = -1/T * deltaH
	     *
	     * Computationally, if Pacc is going to be small, we
	     * can assume Pacc is 0 without actually working it
	     * out.
	     * i.e.,
	     *     if ln Pacc < ln eps, let Pacc = 0
	     * substituting,
	     *     if -deltaH / T < ln eps, let Pacc = 0
	     * rearranging,
	     *     if -deltaH < T * ln eps, let Pacc = 0
	     * This lets us work out whether Pacc will be very
	     * close to zero without dividing anything by T. This
	     * should prevent overflow. Since T is no less
	     * than eps and ln eps is going to have greater
	     * magnitude than eps, underflow when calculating
	     * T * ln eps is not possible. */
	    if (-deltah < t * log_wrapper(LVB_EPS))
	    {
		pacc = 0.0;
		/* Call uni() even though its not required. It
		 * would have been called in LVB 1.0A, so this
		 * helps make results identical to results with
		 * that version. */
		(void) uni();
	    }
	    else	/* possibly accept the change */
	    {
		pacc = exp_wrapper(-deltah/t);
		if (uni() < pacc)	/* do accept the change */
		{
		    probaccd = LVB_TRUE;
		    treeswap(&x, &root, &xdash, &rootdash);
		}
	    }
	    if (probaccd == LVB_TRUE)
	    {
		prev_len = len;
		len = lendash;
	    }
	}
	proposed++;
	if (newtree == LVB_TRUE)
	    accepted++;

	/* decide whether to reduce temperature */
	if (accepted >= maxaccept)	/* enough new trees */
	{
	    failedcnt = 0;  /* this temperature a 'success' */
	    dect = LVB_TRUE;
	}
	else if (proposed >= maxpropose)	/* enough proposals */
	{
	    failedcnt++;
	    if (failedcnt >= maxfail)	/* system frozen */
		break;	/* end of cooling */
	    else	/* decrease temp. */
		dect = LVB_TRUE;
	}

	if (dect == LVB_TRUE)
	{
	    t_n++;	/* originally n is 0 */
	    /* near the start of the function we ensure t1_to_t0 is
	     * either zero or no less than LVB_EPS */
	    if (t1_to_t0 < LVB_EPS)
		t = LVB_EPS;
	    else
	    {
		ln_t = ((double) t_n) * log_wrapper(t1_to_t0) + log_wrapper(t0);
		if (ln_t < log_wrapper(LVB_EPS))
		    t = LVB_EPS;
		else
		    t = pow_wrapper(t1_to_t0, (double) t_n) * t0;
	    }
	    proposed = 0;
	    accepted = 0;
	    dect = LVB_FALSE;
	}
	iter++;
    }

    /* free "local" dynamic heap memory */
    free(x);
    free(xdash);

    return lenbest;

} /* end anneal() */
