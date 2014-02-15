/* LVB (c) Copyright 2003-2006 by Daniel Barker.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** lvb.h - main header for lvb ********** */

/* $Id: lvb.h,v 1.94 2006/02/06 20:09:31 db60 Exp $ */

#ifndef LVB_LVB_H
#define LVB_LVB_H

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "myuni.h"
#include "mymaths.h"
#include "selfconf.h"

/* the program */
#define PROGNAM "lvb"				/* program file name */
#define LVB_VERSION "2.2"			/* version of program */
#define LVB_SUBVERSION "(7th February 2006)"	/* version details e.g. date */

/* DNA bases: bits to set in statesets */
#define A_BIT (1U << 0)
#define C_BIT (1U << 1)
#define G_BIT (1U << 2)
#define T_BIT (1U << 3)
#define O_BIT (1U << 4)

/* values some people may feel the dangerous urge to change */
#define LVB_FNAMSIZE 2000		/* maximum bytes for file names */
#define LVB_INPUTSTRING_SIZE 2000	/* max. bytes for interactive input */
#define UNSET (-1)			/* value of integral vars when unset */
#define STAT_LOG_INTERVAL 100000	/* min. interval for progress log */
#define REROOT_INTERVAL 50		/* change root every ... updates */

/* limits that could be changed but, if increased enormously, might lead to
 * some trouble at some point */
#define MAX_N 100000	/* max. rows */
#define MAX_M 5000000	/* max. cols */

/* implementation-independent limits */
#define LVB_EPS 1E-8		/* 0.0 < DBL_EPSILON < LVB_EPS */
#define MIN_M 1L		/* min. no. of characters for any analysis */
#define MAX_BRANCHES (2 * MAX_N - 3)	/* max. branches per tree */
#define MIN_BRANCHES (2 * MIN_N - 3)	/* max. branches per tree */
#define MIN_N 5L		/* min. no. of objs, for rearrangeable tree */
#define MAX_ALLOC ((size_t) (INT_MAX - 2))	/* max. bytes per dyn. alloc. */
#define MAXSTATES 5		/* max. "true" states in data matrix */

/* limit that could be changed but is likely to be enough */
#define MAX_BOOTSTRAPS 1000000	/* max. bootstrap replicates */

/* unchangeable types */
typedef enum { LVB_FALSE, LVB_TRUE } Lvb_bool;	/* boolean type */

/* matrix and associated information */
typedef struct data
{
    char **row;		/* array of row strings */
    long m;		/* number of columns */
    long n;		/* number of rows */
    char **rowtitle;	/* array of row title strings */ 
} *Dataptr;

/* branch of tree */
typedef struct
{
    long parent;		/* parent branch number, UNSET in root */
    long left;			/* index of first child in tree array */
    long right;			/* index of second child in tree array */
    long object;		/* object number if leaf, otherwise UNSET */
    long changes;		/* changes associated with this branch */
    unsigned char *sset;	/* statesets for all sites */

} Branch;

/* tree stacks */
typedef struct
{
    Branch *tree;	/* pointer to first branch in tree array */
    long root;		/* root of tree */
} Treestack_element;
typedef struct
{
    long size;			/* number of trees currently allocated for */
    long next;			/* next unused element of stack */
    Treestack_element *stack;	/* pointer to first element in stack */
} Treestack;

/* user-, programmer- or program-configurable parameters */
typedef struct
{
    int seed;			/* seed for random number generator */
    long verbose;		/* verboseness level */
    long runs;			/* number of runs */
    long trees;			/* approximate updates to consider */
    long bootstraps;		/* number of bootstrap replicates */
    double t0;			/* initial temperature */
    double t1;			/* second temperature */
    long maxaccept;		/* max. changes at one temperature */
    long maxpropose;		/* max. changes considered at a temperature */
    long maxfail;		/* max. consecutive temps w/ no change */
    Lvb_bool interleaved;	/* LVB_TRUE if matrix is interleaved */
    Lvb_bool fifthstate;	/* if LVB_TRUE, '-' is 'O'; otherwise is '?' */
} Params;

/* fixed file names */
#define MATFNAM "infile"	/* matrix file name */
#define OUTTREEFNAM "outtree"	/* overall best trees */

/* verbose-mode file name bases (run-specific suffixes will be used) */
#define LENFNAM "stat"		/* current tree and length file name prefix */
#define RESFNAM "res"		/* results file name prefix */
#define SUMFNAM "sum"		/* summary of trees per run file name */
#define TREE1FNAM "ini"		/* initial tree file name prefix */

/* default self-configured S.A., hillclimb and randomwalk parameter */
#define TREES_FAST 30000
#define TREES_SLOW 500000

/* assert-like macro, differing in that it writes to standard output,
 * calls crash() not abort(), and works whether or not NDEBUG is defined */
#define lvb_assert(test) \
 ((void) ((test) || (lvb_assertion_fail(#test, __FILE__, __LINE__), 0)))

/* PHYLIP global data */
extern long chars;	/* defined in dnapars.c */

/* LVB global functions */
void *alloc(const size_t bytes, const char *const msg);
long anneal(Treestack *, const Branch *const, long, const double, const double,
 const long, const long, const long, FILE *const, unsigned char **, long, long,
 const long *, long *, Lvb_bool);
long arbreroot(Branch *const, const long);
long brcnt(long);
long childadd(Branch *const, const long, const long);
long cistrcmp(const char *const, const char *const);
Lvb_bool cleanup(void);
void clnclose(FILE *const fp, const char *const fnam);
FILE *clnopen(const char *const nam, const char *const mod);
void clnremove(const char *const fnam);
void crash(const char *const, ...);
void dna_makebin(const Dataptr, Lvb_bool, unsigned char **);
void dnapars_wrapper(void);
char *f2str(FILE *const stream);
Lvb_bool file_exists(const char *const);
void get_bootstrap_weights(long *, long, long);
long getminlen(const Dataptr);
void getparam(Params *);
long getplen(Branch *, const long, const long, const long, const long *);
double get_predicted_length(double, double, long, long, long,
 long);
double get_predicted_trees(double, double, long, long, long,
 long);
long getroot(const Branch *const);
void lvb_assertion_fail(const char *, const char *, int);
void lvb_initialize(void);
Dataptr lvb_matrin(const char *);
void lvb_treeprint (FILE *const stream, const Branch *const barray,
 const long root);
Dataptr matalloc(const long n);
void matchange(Dataptr matrix, const Params rcstruct, const Lvb_bool verbose);
Dataptr matrin(const char *const fnam);
void mutate_spr(Branch *const, const Branch *const, long root);
void mutate_nni(Branch *const, const Branch *const, long);
char *nextnonwspc(const char *string);
void nodeclear(Branch *const, const long);
long objreroot(Branch *const, const long, const long);
void params_change(Params *);
Dataptr phylip_dna_matrin(Lvb_bool);
void phylip_mat_dims_in(long *, long *);
void randtree(Branch *const);
long randpint(const long);
void rowfree(Dataptr matrix);
char *salloc(const long, const char *const);
void scream(const char *const, ...);
void ss_init(Branch *, unsigned char **, long, long);
char *supper(char *const s);
Branch *treealloc(long, long);
void treeclear(Branch *const);
void treecopy(Branch *const, const Branch *const);
long treecmp(const Branch *const, const long, const Branch *const, long);
void treedump(FILE *const, const Branch *const);
void treestack_clear(Treestack *);
long treestack_cnt(Treestack);
long treestack_dump(Treestack *sp, FILE *const outfp);
void treestack_free(Treestack *);
Treestack treestack_new(void);
void treestack_transfer(Treestack *, Treestack *);
long treestack_pop(Branch *, long *, Treestack *);
long treestack_print(Treestack *, FILE *const);
long treestack_push(Treestack *, const Branch *const, const long);
void treeswap(Branch **const, long *const, Branch **const, long *const);

#endif /* LVB_LVB_H */
