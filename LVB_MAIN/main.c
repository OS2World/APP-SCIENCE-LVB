/* LVB (c) Copyright 2003-2006 by Daniel Barker.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** main.c - LVB ********** */

#include "lvb.h"

static const char *rcsid
 = "$Id: main.c,v 1.73 2006/02/06 19:55:46 db60 Exp $";

static Treestack bstack_overall;	/* overall best tree stack */

static void check_stdout(void)
/* Flush standard output, and crash verbosely on error. */
{
    if (fflush(stdout) == EOF)
        crash("write error on standard output");	/* may not work! */
    if (ferror(stdout))
	crash("file error on standard output");		/* may not work! */
}	/* end check_stdout() */

static void smessg(long run)
/* print run start message */
{
    printf("Beginning search %ld\n", run);
    check_stdout();

} /* end smessg() */

static void writeinf(Params prms)
/* write initial details to standard output */
{
    long estimated_trees;	/* estimated trees examined in annealing */
    long estimated_length;	/* estimated best length after annealing */

    printf("\n");

    /* configurable details */
    printf("format               = ");
    if (prms.interleaved == LVB_TRUE)
	printf("INTERLEAVED\n");
    else if (prms.interleaved == LVB_FALSE)
	printf("SEQUENTIAL\n");
    else
	lvb_assert(0);
    printf("treatment of '-'     = ");
    if (prms.fifthstate == LVB_TRUE)
	printf("FIFTH STATE\n");
    else
	printf("UNKNOWN\n");
    printf("seed                 = %d\n", prms.seed);
    printf("duration             = ");
    if (prms.trees == TREES_FAST)
	printf("FAST\n");
    else
	printf("SLOW\n");
    printf("bootstrap replicates = %ld\n", prms.bootstraps);

    if (prms.verbose == LVB_TRUE) {
	printf("trees      = %ld\n", prms.trees);
	printf("verbose    = %d\n", (int) prms.verbose);
	printf("t0         = %g\n", prms.t0);
	printf("t1         = %g\n", prms.t1);
	printf("maxaccept  = %ld\n", prms.maxaccept);
	printf("maxpropose = %ld\n", prms.maxpropose);
	printf("maxfail    = %ld\n", prms.maxfail);
	printf("runs       = %ld\n", prms.runs);

	/* LVB's deductions-in-advance for simulated annealing */
	/* get estimated trees examined and length of best tree,
	 * rounded to the nearest integer */
	estimated_trees = (long) (get_predicted_trees(prms.t0, prms.t1,
	 prms.maxpropose, prms.maxaccept, prms.maxfail, prms.runs) + 0.5);
	 estimated_length = (long) (get_predicted_length(prms.t0,
	 prms.t1, prms.maxpropose, prms.maxaccept,
	     prms.maxfail, prms.runs) + 0.5);

	printf("estimated_trees = %ld\n", estimated_trees);
	printf("estimated_length = %ld\n", estimated_length);

	printf("\n");
    }

} /* end writeinf() */

static void logtree1(const Branch *const barray, const long run, long root)
/* log initial tree for run run (in barray) to outfp */
{
    static char outfnam[LVB_FNAMSIZE]; 	/* current file name */
    int fnamlen;			/* length of current file name */
    FILE *outfp;			/* output file */

    fnamlen = sprintf(outfnam, "%s%ld", TREE1FNAM, run);
    lvb_assert(fnamlen < LVB_FNAMSIZE);	/* shut door if horse bolted */

    /* create tree file */
    outfp = clnopen(outfnam, "w");
    lvb_treeprint(outfp, barray, root);
    clnclose(outfp, outfnam);

} /* end logtree1() */

static long getsoln(Params rcstruct, const long *weight_arr, long *iter_p,
 Lvb_bool log_progress)
/* get and output solution(s) according to parameters in rcstruct;
 * return length of shortest tree(s) found, using weights in weight_arr */
{
    extern Dataptr matrix;		/* data matrix */
    static char fnam[LVB_FNAMSIZE];	/* current file name */
    int fnamlen;			/* length of current file name */
    int i;				/* loop counter */
    long run;				/* current run number */
    long treec;				/* number of trees found */
    long treelength = -1;		/* length of each tree found */
    long best_treelength = LONG_MAX;	/* length of shortest tree found */
    long initroot;			/* initial tree's root */
    FILE *sumfp;			/* best length file */
    FILE *resfp;			/* results file */
    Branch *tree;			/* initial tree */
    Branch *user_tree_ptr = NULL;	/* user-specified initial tree */
    static unsigned char *enc_mat[MAX_N] = { NULL };	/* encoded data mat. */
    Treestack bstack_current;		/* current run's best tree stack */

    /* dynamic "local" heap memory */
    tree = treealloc(brcnt(matrix->n), matrix->m);

    /* Allocation of the initial encoded matrix is non-contiguous because
     * this matrix isn't used much, so any performance penalty won't matter. */
    for (i = 0; i < matrix->n; i++)
        enc_mat[i] = alloc(matrix->m * sizeof(unsigned char), "state sets");
    
    dna_makebin(matrix, rcstruct.fifthstate, enc_mat);

    /* open and entitle statistics file shared by all runs */
    if (rcstruct.verbose == LVB_TRUE) {
	sumfp = clnopen(SUMFNAM, "w");
	fprintf(sumfp, "Run\tInitLen\tBestLen\tTrees\n");
    }
    else
    {
        sumfp = NULL;
    }

    for (run = 0; run < rcstruct.runs; run++)
    {
        if ((run == 0) || (run < (rcstruct.runs - 1)))	/* do from scratch */
	{
	    randtree(tree);
	    ss_init(tree, enc_mat, brcnt(matrix->n), matrix->m);
	    initroot = 0;
	}
	else	/* copy a tree from the overall best stack and start there */
	{
	    /* we achieve "copy" by means of pop followed by re-push */
	    treestack_pop(tree, &initroot, &bstack_overall);
	    treestack_push(&bstack_overall, tree, initroot);
	}

	/* "block local" dynamic heap memory */
        bstack_current = treestack_new();

        if (rcstruct.verbose)
	    smessg(run);
	check_stdout();

	/* start run's entry in sum file */
        if(rcstruct.verbose == LVB_TRUE)
        {
	    fprintf(sumfp, "%ld\t%ld\t", run, getplen(tree, initroot,
	     matrix->m, matrix->n, weight_arr));
	    logtree1(tree, run, initroot);
        }

	/* find solution(s) */
	treelength = anneal(&bstack_current, tree, initroot, rcstruct.t0,
	 rcstruct.t1, rcstruct.maxaccept, rcstruct.maxpropose,
	 rcstruct.maxfail, stdout, enc_mat, matrix->m, matrix->n, weight_arr,
	 iter_p, log_progress);

	/* log this run's solution and its details */
        if (rcstruct.verbose == LVB_TRUE)
        {
	    fnamlen = sprintf(fnam, "%s%ld", RESFNAM, run);
	    lvb_assert(fnamlen < LVB_FNAMSIZE);	/* really too late */
	    resfp = clnopen(fnam, "w");
	    treec = treestack_print(&bstack_current, resfp);
	    clnclose(resfp, fnam);
	    fprintf(sumfp, "%ld\t%ld\n", treelength, treec);

	    /* won't use length summary file until end of next run */
	    fflush(sumfp);
	    if (ferror(sumfp))
	    {
		crash("write error on file %s", SUMFNAM);
	    }
        }

	/* store this run's solution if it's good enough */
        if (treelength < best_treelength)
        {
            treestack_clear(&bstack_overall);
            treestack_transfer(&bstack_overall, &bstack_current);
        }
        else if (treelength == best_treelength)
        {
            treestack_transfer(&bstack_overall, &bstack_current);
        }

        if (rcstruct.verbose == LVB_TRUE)
	    printf("Ending search %ld\n", run);
	check_stdout();

        /* "block-local" dynamic heap memory */
        treestack_free(&bstack_current);
    }
    if (rcstruct.verbose == LVB_TRUE)
	clnclose(sumfp, SUMFNAM);

    /* "local" dynamic heap memory */
    free(tree);
    if (user_tree_ptr != NULL)
	free(user_tree_ptr);

    return treelength;

} /* end getsoln() */

static void logstim(void)
/* log start time with message */
{
    time_t tim;	/* time */

    tim = time(NULL);
    printf("Starting at: %s", ctime(&tim));

} /* end logstim() */

int main(void)
{
    extern Dataptr matrix;	/* data matrix */
    int val;			/* return value */
    Params rcstruct;		/* configurable parameters */
    long i;			/* loop counter */
    long m;			/* sites per sequence */
    long n;			/* sequences in the data matrix */
    long iter;			/* iterations of annealing algorithm */
    long replicate_no = 0;	/* current bootstrap replicate number */
    long trees_found_total = 0;	/* number of trees found, overall */
    long trees_found;		/* number of trees for current replicate */
    double total_iter = 0.0;	/* total iterations across all replicates */
    long final_length;		/* length of shortest tree(s) found */
    long m_including_constcols;	/* site count before constant sites removed */
    FILE *outtreefp;		/* best trees found overall */
    static long weight_arr[MAX_M];	/* weights for sites */
    Lvb_bool log_progress;	/* whether or not to log anneal search */

    /* global files */

    /* entitle standard output */
    printf("\n* This is %s version %s %s, written by Daniel Barker *\n",
     PROGNAM, LVB_VERSION, LVB_SUBVERSION);
    printf("\n");
    printf("Reference:\n");
    printf("Barker, D. 2004. LVB: Parsimony and simulated annealing in the "
    "search\n"
    "for phylogenetic trees. Bioinformatics, 20, 274-275.\n");
    printf("\n");

    lvb_initialize();

    getparam(&rcstruct);
    phylip_mat_dims_in(&n, &m);
    logstim();
    params_change(&rcstruct);

    matrix = phylip_dna_matrin(rcstruct.interleaved);
    lvb_assert((matrix->m == m) && (matrix->n == n));

    /* "file-local" dynamic heap memory: set up best tree stacks */
    bstack_overall = treestack_new();

    writeinf(rcstruct);
    m_including_constcols = matrix->m;
    matchange(matrix, rcstruct, rcstruct.verbose);	/* cut columns */

    if (rcstruct.verbose == LVB_TRUE)
    {
	printf("getminlen: %ld\n", getminlen(matrix));
    }
    rinit(rcstruct.seed);
    if (rcstruct.bootstraps > 0)
    {
	log_progress = LVB_FALSE;
	printf("\nReplicate:      Rearrangements: Trees:          Length:\n");
    }
    else
	log_progress = LVB_TRUE;
    outtreefp = clnopen(OUTTREEFNAM, "w");
    do
    {
	iter = 0;
	if (rcstruct.bootstraps > 0)
	{
	    get_bootstrap_weights(weight_arr, matrix->m,
		m_including_constcols - matrix->m);
	}
	else
	{
	    for (i = 0; i < matrix->m; i++)
		weight_arr[i] = 1;
	}

	final_length = getsoln(rcstruct, weight_arr, &iter, log_progress);
	trees_found = treestack_print(&bstack_overall, outtreefp);
	trees_found_total += trees_found;
	treestack_clear(&bstack_overall);
	replicate_no++;
	if (rcstruct.bootstraps > 0)
	{
	    printf("%-16ld%-16ld%-16ld%ld\n", replicate_no, iter, trees_found,
		final_length);
	    total_iter += (double) iter;
	} else
	    printf("\nRearrangements tried: %ld\n", iter);
    } while (replicate_no < rcstruct.bootstraps);
    clnclose(outtreefp, OUTTREEFNAM);

    printf("\n");

    if (rcstruct.bootstraps > 0)
	printf("Total rearrangements tried across all replicates: %g\n\n",
	 total_iter);

    if (trees_found_total == 1L)
    {
	printf("1 most parsimonious tree of length %ld written to file "
	    "'%s'\n", final_length, OUTTREEFNAM);
    }
    else
    {
	if (rcstruct.bootstraps > 0)
	    printf("%ld trees written to file '%s'\n", trees_found_total,
		OUTTREEFNAM);
	else
	{
	    printf("%ld equally parsimonious trees of length %ld written to "
         "file '%s'\n", trees_found_total, final_length, OUTTREEFNAM);
	}
    }

    rowfree(matrix);

    /* "file-local" dynamic heap memory */
    treestack_free(&bstack_overall);

    if (cleanup() == LVB_TRUE)
	val = EXIT_FAILURE;
    else
	val = EXIT_SUCCESS;

    return val;

} /* end main() */
