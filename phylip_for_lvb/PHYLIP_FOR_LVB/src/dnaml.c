
#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1993-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   Dan Fineman, and Patrick Colacurcio.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

typedef struct valrec {
  double rat, ratxi, ratxv, orig_zz, z1, y1, z1zz, z1yy, xiz1,
	 xiy1xv;
  double *ww, *zz, *wwzz, *vvzz; 
} valrec;

typedef long vall[maxcategs];
typedef double contribarr[maxcategs];

#ifndef OLDC
/* function prototypes */
void   dnamlcopy(tree *, tree *, short, short);
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   inputoptions(void);
void   makeweights(void);
void   getinput(void);
void   inittable_for_usertree(FILE *, char *, char *);
void   inittable(void);
double evaluate(node *, boolean);

void   alloc_nvd (long, nuview_data *);
void   free_nvd (long, nuview_data *);
void   nuview(node *);
void   slopecurv(node *, double, double *, double *, double *);
void   makenewv(node *);
void   update(node *);
void   smooth(node *);
void   insert_(node *, node *, boolean);
void   dnaml_re_move(node **, node **);
void   buildnewtip(long, tree *);

void   buildsimpletree(tree *);
void   addtraverse(node *, node *, boolean);
void   rearrange(node *, node *);
void   initdnamlnode(node **, node **, node *, long, long, long *, long *,
		 initops, pointarray, pointarray, Char *, Char *, FILE *);
void   dnaml_coordinates(node *, double, long *, double *);
void   dnaml_printree(void);
void   sigma(node *, double *, double *, double *);
void   describe(node *);
void   reconstr(node *, long );
void   rectrav(node *, long, long);

void   summarize(void);
void   dnaml_treeout(node *);
void   inittravtree(node *);
void   treevaluate(void);
void   maketree(void);
void   clean_up(void);
/* function prototypes */
#endif


extern FILE *infile, *outfile, *intree, *outtree, *catfile, *weightfile;
extern long spp, outgrno, endsite, which;
extern boolean ibmpc, ansi, interleaved, printdata, outgropt, treeprint;
extern naym *nayme;
extern sequence y;
extern steptr weight, category, alias, location, ally;

double fracchange;
long rcategs;
boolean haslengths;

Char infilename[100], outfilename[100], intreename[100], outtreename[100],
     catfilename[100], weightfilename[100];
double *rate, *rrate, *probcat;
long nonodes2, sites, weightsum, categs,
  datasets, ith, njumble, jumb;
long inseed, parens;
boolean  freqsfrom, global, jumble, weights, trout, usertree, reconsider,
  ctgry, rctgry, auto_, hypstate, ttr, progress, mulsets, justwts,
  firstset, improve, smoothit, polishing, lngths;
tree curtree, bestree, bestree2, priortree;
node *qwhere, *grbg;
double xi, xv, ttratio, ttratio0, freqa, freqc, freqg, freqt, freqr,
  freqy, freqar, freqcy, freqgr, freqty, 
  lambda, bestyet;
long *enterorder;
steptr aliasweight;
contribarr *contribution, like, nulike, clai;
double **term, **slopeterm, **curveterm;
longer seed;
Char progname[40];
char basechar[16]="acmgrsvtwyhkdbn";

/* Local variables for maketree, propagated globally for c version: */
long k, nextsp, numtrees, maxwhich, mx, mx0;
double dummy, maxlogl;
boolean succeeded, smoothed;
double **l0gf;
double *l0gl;
valrec ***tbl;
Char ch, ch2;
long col;
vall *mp;


void dnamlcopy(tree *a, tree *b, short nonodes, short categs)
{
  /* used in dnaml & dnamlk */
  short i, j=0;
  node *p, *q;

  for (i = 0; i < spp; i++) {
    copynode(a->nodep[i], b->nodep[i], categs);
    if (a->nodep[i]->back) {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1];
      else if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1]->next)
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next;
      else
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next->next;
    }
    else b->nodep[i]->back = NULL;
  }
  for (i = spp; i < nonodes; i++) {
    p = a->nodep[i];
    q = b->nodep[i];
    for (j = 1; j <= 3; j++) {
      copynode(p, q, categs);
      if (p->back) {
        if (p->back == a->nodep[p->back->index - 1])
          q->back = b->nodep[p->back->index - 1];
        else if (p->back == a->nodep[p->back->index - 1]->next)
          q->back = b->nodep[p->back->index - 1]->next;
        else
          q->back = b->nodep[p->back->index - 1]->next->next;
      }
      else
        q->back = NULL;
      p = p->next;
      q = q->next;
    }
  }
  b->likelihood = a->likelihood;
  b->start = a->start;               /* start used in dnaml only */
  b->root = a->root;                 /* root used in dnamlk only */
}  /* dnamlcopy plc*/


void getoptions()
{
  /* interactively set options */
  long inseed0;
  Char ch;
  boolean didchangecat, didchangercat;
  double probsum;

  fprintf(outfile, "\nNucleic acid sequence Maximum Likelihood");
  fprintf(outfile, " method, version %s\n\n",VERSION);
  putchar('\n');
  ctgry = false;
  didchangecat = false;
  rctgry = false;
  didchangercat = false;
  categs = 1;
  rcategs = 1;
  auto_ = false;
  freqsfrom = true;
  global = false;
  hypstate = false;
  improve = false;
  jumble = false;
  njumble = 1;
  lngths = false;
  lambda = 1.0;
  outgrno = 1;
  outgropt = false;
  reconsider = false;
  trout = true;
  ttratio = 2.0;
  ttr = false;
  usertree = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  interleaved = true;
  for (;;){
    cleerhome();
    printf("Nucleic acid sequence Maximum Likelihood");
    printf(" method, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
	   (usertree ? "No, use user trees in input file" : "Yes"));
    if (usertree) {
      printf("  L          Use lengths from user trees?  %s\n",
	       (lngths ? "Yes" : "No"));
    }
    printf("  T        Transition/transversion ratio:%8.4f\n",
	   (ttr ? ttratio : 2.0));
    printf("  F       Use empirical base frequencies?  %s\n",
	   (freqsfrom ? "Yes" : "No"));
    printf("  C                One category of sites?");
    if (!ctgry || categs == 1)
      printf("  Yes\n");
    else
      printf("  %ld categories of sites\n", categs);
    printf("  R     One region of substitution rates?");
    if (!rctgry || rcategs == 1)
      printf("  Yes\n");
    else {
      printf("  %ld categories of regions\n", rcategs);
      printf("  A   Rates at adjacent sites correlated?");
      if (!auto_)
	printf("  No, they are independent\n");
      else
	printf("  Yes, mean block length =%6.1f\n", 1.0 / lambda);
    }
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    if ((!usertree) || reconsider) {
      printf("  S        Speedier but rougher analysis?  %s\n",
             (improve ? "No, not rough" : "Yes"));
      printf("  G                Global rearrangements?  %s\n",
	     (global ? "Yes" : "No"));
    }
    if (!usertree) {
      printf("  J   Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
	printf("  No. Use input order\n");
    }
    else
       printf("  V    Rearrange starting with user tree?  %s\n",
             (reconsider ? "Yes" : "No"));
    printf("  O                        Outgroup root?  %s%3ld\n",
	   (outgropt ? "Yes, at sequence number" :
                       "No, use as outgroup species"),outgrno);
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", datasets,
               (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
	   (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
	   (ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
	   (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
	   (treeprint ? "Yes" : "No"));
    printf("  4       Write out trees onto tree file?  %s\n",
	   (trout ? "Yes" : "No"));
    printf("  5   Reconstruct hypothetical sequences?  %s\n",
           (hypstate ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("ULTFCRAWSGJVOMI012345",ch) != NULL){
      switch (ch) {

      case 'F':
	freqsfrom = !freqsfrom;
	if (!freqsfrom) {
	  initfreqs(&freqa, &freqc, &freqg, &freqt);
	}
	break;
	
      case 'C':
	ctgry = !ctgry;
	if (ctgry) {
          printf("\nSitewise user-assigned categories:\n\n");
          initcatn(&categs);
          if (rate){
	    free(rate);
	  }
          rate    = (double *) calloc (1,categs * sizeof(double));
          didchangecat = true;
          initcategs(categs, rate);
        }
        break;

      case 'R':
	rctgry = !rctgry;
	if (!rctgry)
	  auto_ = false;
	if (rctgry) {
          printf("\nRegional rates:\n\n");
          initcatn(&rcategs);
          if (probcat){
            free(probcat);
	    free(rrate);
	  }
          probcat = (double *) calloc (1,rcategs * sizeof(double));
          rrate   = (double *) calloc (1,rcategs * sizeof(double));
          didchangercat = true;
          initcategs(rcategs, rrate);
          initprobcat(rcategs, &probsum, probcat);
	}
	break;
	
      case 'A':
	auto_ = !auto_;
	if (auto_)
	  initlambda(&lambda);
	break;
	
      case 'W':
        weights = !weights;
        break;

      case 'S':
        improve = !improve;
        break;

      case 'G':
	global = !global;
	break;
	
      case 'J':
	jumble = !jumble;
	if (jumble)
	  initjumble(&inseed, &inseed0, seed, &njumble);
	else njumble = 1;
	break;
	
      case 'L':
	lngths = !lngths;
	break;
	
      case 'O':
	outgropt = !outgropt;
	if (outgropt)
	  initoutgroup(&outgrno, spp);
	break;
	
      case 'T':
	ttr = !ttr;
	if (ttr) {
	  initratio(&ttratio);
	}
	break;
	
      case 'U':
	usertree = !usertree;
	break;

      case 'V':
        reconsider = !reconsider;
        break;

      case 'M':
	mulsets = !mulsets;
	if (mulsets) {
          printf("Multiple data sets or multiple weights?");
          do {
            printf(" (type D or W)\n");
            scanf("%c%*[^\n]", &ch2);
            uppercase(&ch2);
          } while ((ch2 != 'W') && (ch2 != 'D'));
          justwts = (ch2 == 'W');
          if (justwts)
            justweights(&datasets);
          else
	    initdatasets(&datasets);
          if (!jumble) {
            jumble = true;
	    initjumble(&inseed, &inseed0, seed, &njumble);
          }
        }
	break;
	
      case 'I':
	interleaved = !interleaved;
	break;
	
      case '0':
	initterminal(&ibmpc, &ansi);
	break;
	
      case '1':
	printdata = !printdata;
	break;
	
      case '2':
	progress = !progress;
	break;
	
      case '3':
	treeprint = !treeprint;
	break;
	
      case '4':
	trout = !trout;
	break;

      case '5':
	hypstate = !hypstate;
	break;
      }
    } else
      printf("Not a possible option!\n");
  }
  if (!didchangercat){
    rrate      = (double *) calloc (1,rcategs*sizeof(double));
    probcat    = (double *) calloc (1,rcategs*sizeof(double));
    rrate[0]   = 1.0;
    probcat[0] = 1.0;
  }
  if (!didchangecat){
    rate       = (double *) calloc (1,categs*sizeof(double));
    rate[0]    = 1.0;
  }
}  /* getoptions */


void allocrest()
{
  long i;

  y = (Char **) calloc (1,spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    y[i] = (Char *) calloc (1,sites*sizeof(Char));
  nayme       = (naym *) calloc (1,spp*sizeof(naym));;
  enterorder  = (long *) calloc (1,spp*sizeof(long));
  category    = (long *) calloc (1,sites*sizeof(long));
  weight      = (long *) calloc (1,sites*sizeof(long));
  alias       = (long *) calloc (1,sites*sizeof(long));
  ally        = (long *) calloc (1,sites*sizeof(long));
  location    = (long *) calloc (1,sites*sizeof(long));
  aliasweight = (long *) calloc (1,sites*sizeof(long));
}  /* allocrest */


void doinit()
{ /* initializes variables */

  inputnumbers(&spp, &sites, &nonodes2, 2);
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, sites);
  alloctree(&curtree.nodep, nonodes2, usertree);
  allocrest();
  if (usertree && !reconsider)
    return;
  alloctree(&bestree.nodep, nonodes2, 0);
  alloctree(&priortree.nodep, nonodes2, 0);
  if (njumble <= 1)
    return;
  alloctree(&bestree2.nodep, nonodes2, 0);
}  /* doinit */


void inputoptions()
{
  long i;

  if (!firstset)
    samenumsp(&sites, ith);
  if (firstset) {
    for (i = 0; i < sites; i++)
      category[i] = 1;
    for (i = 0; i < sites; i++)
      weight[i] = 1;
  }
  if (justwts || weights)
    inputweights(sites, weight, &weights);
  weightsum = 0;
  for (i = 0; i < sites; i++)
    weightsum += weight[i];
  if ((ctgry && categs > 1) && (firstset || !justwts)) {
    inputcategs(0, sites, category, categs, "DnaML");
    if (printdata)
      printcategs(outfile, sites, category, "Site categories");
  }
  if (weights && printdata)
    printweights(outfile, 0, sites, weight, "Sites");
}  /* inputoptions */


void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= sites; i++) {
    alias[i - 1] = i;
    ally[i - 1] = 0;
    aliasweight[i - 1] = weight[i - 1];
    location[i - 1] = 0;
  }
  sitesort2   (sites, aliasweight);
  sitecombine2(sites, aliasweight);
  sitescrunch2(sites, 1, 2, aliasweight);
  for (i = 1; i <= sites; i++) {
    if (aliasweight[i - 1] > 0)
      endsite = i;
  }
  for (i = 1; i <= endsite; i++) {
    location[alias[i - 1] - 1] = i;
    ally[alias[i - 1] - 1] = alias[i - 1];
  }
  term = (double **) calloc (1,endsite * sizeof(double *));
  for (i = 0; i < endsite; i++)
     term[i] = (double *) calloc (1,rcategs * sizeof(double));
  slopeterm = (double **) calloc (1,endsite * sizeof(double *));
  for (i = 0; i < endsite; i++)
     slopeterm[i] = (double *) calloc (1,rcategs * sizeof(double));
  curveterm = (double **) calloc (1,endsite * sizeof(double *));
  for (i = 0; i < endsite; i++)
     curveterm[i] = (double *) calloc (1,rcategs * sizeof(double));
  mp = (vall *) calloc (1,sites*sizeof(vall));
  contribution = (contribarr *) calloc (1,endsite*sizeof(contribarr));
}  /* makeweights */


void getinput()
{
  /* reads the input data */
  if (!justwts || firstset)
    inputoptions();
  if (!freqsfrom)
    getbasefreqs(freqa, freqc, freqg, freqt, &freqr, &freqy, &freqar, &freqcy,
		 &freqgr, &freqty, &ttratio, &xi, &xv, &fracchange,
		 freqsfrom, true);
  if (!justwts || firstset)
    inputdata(sites);
  makeweights();
  setuptree2(curtree);
  if (!usertree || reconsider) {
    setuptree2(bestree);
    setuptree2(priortree);
    if (njumble > 1)
      setuptree2(bestree2);
  }
  allocx(nonodes2, rcategs, curtree.nodep, usertree);
  if (!usertree || reconsider) {
    allocx(nonodes2, rcategs, bestree.nodep, 0);
    allocx(nonodes2, rcategs, priortree.nodep, 0);
    if (njumble > 1)
      allocx(nonodes2, rcategs, bestree2.nodep, 0);
  }
  makevalues2(rcategs, curtree.nodep, endsite, spp, y, alias);  
  if (freqsfrom) {
    empiricalfreqs(&freqa, &freqc, &freqg, &freqt, aliasweight,
                    curtree.nodep);
    getbasefreqs(freqa, freqc, freqg, freqt, &freqr, &freqy, &freqar, &freqcy,
		 &freqgr, &freqty, &ttratio, &xi, &xv, &fracchange,
		 freqsfrom, true);
  }
  if (!justwts || firstset)
    fprintf(outfile, "\nTransition/transversion ratio = %10.6f\n\n", ttratio);
}  /* getinput */


void inittable_for_usertree(FILE *intree, char *progname, char *intreename)
{
  /* If there's a user tree, then the ww/zz/wwzz/vvzz elements need
     to be allocated appropriately. */
  long num_comma;
  long i, j;

  /* First, figure out the largest possible furcation, i.e. the number
     of commas plus one */
  countcomma(&intree, &num_comma);
  num_comma++;
  
  for (i = 0; i < rcategs; i++) {
    for (j = 0; j < categs; j++) {
      /* Free the stuff allocated assuming bifurcations */
      free (tbl[i][j]->ww);
      free (tbl[i][j]->zz);
      free (tbl[i][j]->wwzz);
      free (tbl[i][j]->vvzz);

      /* Then allocate for worst-case multifurcations */
      tbl[i][j]->ww   = (double *) calloc (1, num_comma * sizeof (double));
      tbl[i][j]->zz   = (double *) calloc (1, num_comma * sizeof (double));
      tbl[i][j]->wwzz = (double *) calloc (1, num_comma * sizeof (double));
      tbl[i][j]->vvzz = (double *) calloc (1, num_comma * sizeof (double));
    }
  }
}  /* inittable_for_usertree */

void inittable()
{
  /* Define a lookup table. Precompute values and print them out in tables */
  long i, j;
  double sumrates;
  
  tbl = (valrec ***) calloc (1,rcategs * sizeof(valrec **));
  for (i = 0; i < rcategs; i++) {
    tbl[i] = (valrec **) calloc (1,categs*sizeof(valrec *));
    for (j = 0; j < categs; j++)
      tbl[i][j] = (valrec *) calloc (1,sizeof(valrec));
  }

  for (i = 0; i < rcategs; i++) {
    for (j = 0; j < categs; j++) {
      tbl[i][j]->rat = rrate[i]*rate[j];
      tbl[i][j]->ratxi = tbl[i][j]->rat * xi;
      tbl[i][j]->ratxv = tbl[i][j]->rat * xv;

      /* Allocate assuming bifurcations, will be changed later if
         neccesarry (i.e. there's a user tree) */
      tbl[i][j]->ww   = (double *) calloc (1, 2 * sizeof (double));
      tbl[i][j]->zz   = (double *) calloc (1, 2 * sizeof (double));
      tbl[i][j]->wwzz = (double *) calloc (1, 2 * sizeof (double));
      tbl[i][j]->vvzz = (double *) calloc (1, 2 * sizeof (double));
    }
  }
  sumrates = 0.0;
  for (i = 0; i < endsite; i++) {
    for (j = 0; j < rcategs; j++)
      sumrates += aliasweight[i] * probcat[j] 
	* tbl[j][category[alias[i] - 1] - 1]->rat;
  }
  sumrates /= (double)sites;
  for (i = 0; i < rcategs; i++)
    for (j = 0; j < categs; j++) {
      tbl[i][j]->rat /= sumrates;
      tbl[i][j]->ratxi /= sumrates;
      tbl[i][j]->ratxv /= sumrates;
    }
  if (rcategs > 1) {
    fprintf(outfile, "\nRegion type     Rate of change    Probability\n\n");
    for (i = 0; i < rcategs; i++)
      fprintf(outfile, "%9ld%16.3f%17.3f\n", i+1, rrate[i], probcat[i]);
    putc('\n', outfile);
    if (auto_)
      fprintf(outfile,
     "Expected length of a patch of sites having the same rate = %8.3f\n",
             1/lambda);
    putc('\n', outfile);
  }
  if (categs > 1) {
    fprintf(outfile, "\nSite category   Rate of change\n\n");
    for (i = 0; i < categs; i++)
      fprintf(outfile, "%9ld%16.3f\n", i+1, rate[i]);
  }
  if ((rcategs  > 1) || (categs >> 1))
    fprintf(outfile, "\n\n");
}  /* inittable */


double evaluate(node *p, boolean saveit)
{
  contribarr tterm;
  double sum, sum2, sumc, y, lz, y1, z1zz, z1yy, prod12, prod1, prod2, prod3,
	 sumterm, lterm;
  long i, j, k, lai;
  node *q;
  sitelike x1, x2;

  sum = 0.0;
  q = p->back;
  y = p->v;
  lz = -y;
  for (i = 0; i < rcategs; i++)
    for (j = 0; j < categs; j++) {
    tbl[i][j]->orig_zz = exp(tbl[i][j]->ratxi * lz);
    tbl[i][j]->z1 = exp(tbl[i][j]->ratxv * lz);
    tbl[i][j]->z1zz = tbl[i][j]->z1 * tbl[i][j]->orig_zz;
    tbl[i][j]->z1yy = tbl[i][j]->z1 - tbl[i][j]->z1zz;
  }
  for (i = 0; i < endsite; i++) {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++) {
      if (y > 0.0) {
	y1 = 1.0 - tbl[j][k]->z1;
	z1zz = tbl[j][k]->z1zz;
	z1yy = tbl[j][k]->z1yy;
      } else {
	y1 = 0.0;
	z1zz = 1.0;
	z1yy = 0.0;
      }
      memcpy(x1, p->x[i][j], sizeof(sitelike));
      prod1 = freqa * x1[0] + freqc * x1[(long)C - (long)A] +
	     freqg * x1[(long)G - (long)A] + freqt * x1[(long)T - (long)A];
      memcpy(x2, q->x[i][j], sizeof(sitelike));
      prod2 = freqa * x2[0] + freqc * x2[(long)C - (long)A] +
	     freqg * x2[(long)G - (long)A] + freqt * x2[(long)T - (long)A];
      prod3 = (x1[0] * freqa + x1[(long)G - (long)A] * freqg) *
	      (x2[0] * freqar + x2[(long)G - (long)A] * freqgr) +
	  (x1[(long)C - (long)A] * freqc + x1[(long)T - (long)A] * freqt) *
	 (x2[(long)C - (long)A] * freqcy + x2[(long)T - (long)A] * freqty);
      prod12 = freqa * x1[0] * x2[0] +
	       freqc * x1[(long)C - (long)A] * x2[(long)C - (long)A] +
	       freqg * x1[(long)G - (long)A] * x2[(long)G - (long)A] +
	       freqt * x1[(long)T - (long)A] * x2[(long)T - (long)A];
      tterm[j] = z1zz * prod12 + z1yy * prod3 + y1 * prod1 * prod2;
    }
    sumterm = 0.0;
    for (j = 0; j < rcategs; j++)
      sumterm += probcat[j] * tterm[j];
    lterm = log(sumterm);
    for (j = 0; j < rcategs; j++)
      clai[j] = tterm[j] / sumterm;
    memcpy(contribution[i], clai, rcategs*sizeof(double));
    if (saveit && !auto_ && usertree)
      l0gf[which - 1][i] = lterm;
    sum += aliasweight[i] * lterm;
  }
  for (j = 0; j < rcategs; j++)
    like[j] = 1.0;
  for (i = 0; i < sites; i++) {
    sumc = 0.0;
    for (k = 0; k < rcategs; k++)
      sumc += probcat[k] * like[k];
    sumc *= lambda;
    if ((ally[i] > 0) && (location[ally[i]-1] > 0)) {
      lai = location[ally[i] - 1];
      memcpy(clai, contribution[lai - 1], rcategs*sizeof(double));
      for (j = 0; j < rcategs; j++)
        nulike[j] = ((1.0 - lambda) * like[j] + sumc) * clai[j];
    } else {
      for (j = 0; j < rcategs; j++)
        nulike[j] = ((1.0 - lambda) * like[j] + sumc);
    }
    memcpy(like, nulike, rcategs*sizeof(double));
  }
  sum2 = 0.0;
  for (i = 0; i < rcategs; i++)
    sum2 += probcat[i] * like[i];
  sum += log(sum2);
  curtree.likelihood = sum;
  if (!saveit || auto_ || !usertree)
    return sum;
  l0gl[which - 1] = sum;
  if (which == 1) {
    maxwhich = 1;
    maxlogl = sum;
    return sum;
  }
  if (sum > maxlogl) {
    maxwhich = which;
    maxlogl = sum;
  }
  return sum;
}  /* evaluate */


void alloc_nvd (long num_sibs, nuview_data *local_nvd)
{
  /* Allocate blocks of memory appropriate for the number of siblings
     a given node has */
  local_nvd->yy     = (double *) calloc (1, num_sibs * sizeof (double));
  local_nvd->wwzz   = (double *) calloc (1, num_sibs * sizeof (double));
  local_nvd->vvzz   = (double *) calloc (1, num_sibs * sizeof (double));
  local_nvd->vzsumr = (double *) calloc (1, num_sibs * sizeof (double));
  local_nvd->vzsumy = (double *) calloc (1, num_sibs * sizeof (double));
  local_nvd->sum    = (double *) calloc (1, num_sibs * sizeof (double));
  local_nvd->sumr   = (double *) calloc (1, num_sibs * sizeof (double));
  local_nvd->sumy   = (double *) calloc (1, num_sibs * sizeof (double));
  local_nvd->xx     = (sitelike *) calloc (1, num_sibs * sizeof (sitelike));
}  /* alloc_nvd */


void free_nvd (long num_sibs, nuview_data *local_nvd)
{
  /* The natural complement to the alloc version */
  free (local_nvd->yy);
  free (local_nvd->wwzz);
  free (local_nvd->vvzz);
  free (local_nvd->vzsumr);
  free (local_nvd->vzsumy);
  free (local_nvd->sum);
  free (local_nvd->sumr);
  free (local_nvd->sumy);
  free (local_nvd->xx);
}  /* free_nvd */


void nuview(node *p)
{
  long i, j, k, num_sibs, sib_index;
  nuview_data *local_nvd;
  node *sib_ptr, *sib_back_ptr;
  sitelike p_xx;
  double lw;

  /* Figure out how many siblings the current node has */
  num_sibs    = count_sibs (p);

  /* Recursive calls, should be called for all children */
  sib_ptr = p;
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    if (!sib_back_ptr->tip &&
	!sib_back_ptr->initialized)
      nuview (sib_back_ptr);
  }

  /* Allocate the structure and blocks therein for variables used in
     this function */
  local_nvd = (nuview_data *) calloc (1, sizeof (nuview_data));
  alloc_nvd (num_sibs, local_nvd);


  /* Loop 1: makes assignments to tbl based on some combination of
     what's already in tbl and the children's value of v */
  sib_ptr = p;
  for (sib_index=0; sib_index < num_sibs; sib_index++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    
    lw = - (sib_back_ptr->v);

    for (i = 0; i < rcategs; i++)
      for (j = 0; j < categs; j++) {
	tbl[i][j]->ww[sib_index]   = exp(tbl[i][j]->ratxi * lw);
	tbl[i][j]->zz[sib_index]   = exp(tbl[i][j]->ratxv * lw);
	tbl[i][j]->wwzz[sib_index] = tbl[i][j]->ww[sib_index] * tbl[i][j]->zz[sib_index];
	tbl[i][j]->vvzz[sib_index] = (1.0 - tbl[i][j]->ww[sib_index]) *
	  tbl[i][j]->zz[sib_index];
      }
  }
  
  /* Loop 2: */
  for (i = 0; i < endsite; i++) {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++) {

      /* Loop 2.1 */
      sib_ptr = p;
      for (sib_index=0; sib_index < num_sibs; sib_index++) {
	sib_ptr         = sib_ptr->next;
	sib_back_ptr    = sib_ptr->back;

	local_nvd->wwzz[sib_index] = tbl[j][k]->wwzz[sib_index];
	local_nvd->vvzz[sib_index] = tbl[j][k]->vvzz[sib_index];
	local_nvd->yy[sib_index]   = 1.0 - tbl[j][k]->zz[sib_index];
	memcpy(local_nvd->xx[sib_index],
	       sib_back_ptr->x[i][j],
	       sizeof(sitelike));
      }

      /* Loop 2.2 */
      sib_ptr = p;
      for (sib_index=0; sib_index < num_sibs; sib_index++) {
	sib_ptr         = sib_ptr->next;
	sib_back_ptr    = sib_ptr->back;

	local_nvd->sum[sib_index] =
	  local_nvd->yy[sib_index] *
	  (freqa * local_nvd->xx[sib_index][(long)A] +
	   freqc * local_nvd->xx[sib_index][(long)C] +
	   freqg * local_nvd->xx[sib_index][(long)G] +
	   freqt * local_nvd->xx[sib_index][(long)T]);
	local_nvd->sumr[sib_index] =
	  freqar * local_nvd->xx[sib_index][(long)A] +
	  freqgr * local_nvd->xx[sib_index][(long)G];
	local_nvd->sumy[sib_index] =
	  freqcy * local_nvd->xx[sib_index][(long)C] +
	  freqty * local_nvd->xx[sib_index][(long)T];
	local_nvd->vzsumr[sib_index] =
	  local_nvd->vvzz[sib_index] * local_nvd->sumr[sib_index];
	local_nvd->vzsumy[sib_index] =
	  local_nvd->vvzz[sib_index] * local_nvd->sumy[sib_index];
      }

      /* Initialize to one, multiply incremental values for every
         sibling a node has */
      p_xx[(long)A] = 1 ;
      p_xx[(long)C] = 1 ; 
      p_xx[(long)G] = 1 ;
      p_xx[(long)T] = 1 ;

      for (sib_index=0; sib_index < num_sibs; sib_index++) {
	p_xx[(long)A] *=
	  local_nvd->sum[sib_index] +
	  local_nvd->wwzz[sib_index] *
	  local_nvd->xx[sib_index][(long)A] +
	  local_nvd->vzsumr[sib_index];
	p_xx[(long)C] *=
	  local_nvd->sum[sib_index] +
	  local_nvd->wwzz[sib_index] *
	  local_nvd->xx[sib_index][(long)C] +
	  local_nvd->vzsumy[sib_index];
	p_xx[(long)G] *=
	  local_nvd->sum[sib_index] +
	  local_nvd->wwzz[sib_index] *
	  local_nvd->xx[sib_index][(long)G] +
	  local_nvd->vzsumr[sib_index];
	p_xx[(long)T] *=
	  local_nvd->sum[sib_index] +
	  local_nvd->wwzz[sib_index] *
	  local_nvd->xx[sib_index][(long)T] +
	  local_nvd->vzsumy[sib_index];
      }

      /* And the final point of this whole function: */
      memcpy(p->x[i][j], p_xx, sizeof(sitelike));
    }
  }

  p->initialized = true;

  free_nvd (num_sibs, local_nvd);
  free (local_nvd);
}  /* nuview */


void slopecurv(node *p,double y,double *like,double *slope,double *curve)
{
   /* compute log likelihood, slope and curvature at node p */
  long i, j, k, lai;
  double sum, sumc, sumterm, lterm, sumcs, sumcc, sum2, slope2, curve2,
        temp;
  double lz, zz, z1, zzs, z1s, zzc, z1c, aa, bb, cc,
         prod1, prod2, prod12, prod3;
  contribarr thelike, nulike, nuslope, nucurve,
    theslope, thecurve, clai, cslai, cclai;
  node *q;
  sitelike x1, x2;

  q = p->back;
  sum = 0.0;
  lz = -y;
  for (i = 0; i < rcategs; i++)
    for (j = 0; j < categs; j++) {
      tbl[i][j]->orig_zz = exp(tbl[i][j]->rat * lz);
      tbl[i][j]->z1 = exp(tbl[i][j]->ratxv * lz);
    }
  for (i = 0; i < endsite; i++) {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++) {
      if (y > 0.0) {
	zz = tbl[j][k]->orig_zz;
	z1 = tbl[j][k]->z1;
      } else {
	zz = 1.0;
	z1 = 1.0;
      }
      zzs = -tbl[j][k]->rat * zz ;
      z1s = -tbl[j][k]->ratxv * z1 ;
      temp = tbl[j][k]->rat;
      zzc = temp * temp * zz;
      temp = tbl[j][k]->ratxv;
      z1c = temp * temp * z1;
      memcpy(x1, p->x[i][j], sizeof(sitelike));
      prod1 = freqa * x1[0] + freqc * x1[(long)C - (long)A] +
	    freqg * x1[(long)G - (long)A] + freqt * x1[(long)T - (long)A];
      memcpy(x2, q->x[i][j], sizeof(sitelike));
      prod2 = freqa * x2[0] + freqc * x2[(long)C - (long)A] +
	    freqg * x2[(long)G - (long)A] + freqt * x2[(long)T - (long)A];
      prod3 = (x1[0] * freqa + x1[(long)G - (long)A] * freqg) *
	      (x2[0] * freqar + x2[(long)G - (long)A] * freqgr) +
	(x1[(long)C - (long)A] * freqc + x1[(long)T - (long)A] * freqt) *
	(x2[(long)C - (long)A] * freqcy + x2[(long)T - (long)A] * freqty);
      prod12 = freqa * x1[0] * x2[0] +
	       freqc * x1[(long)C - (long)A] * x2[(long)C - (long)A] +
	       freqg * x1[(long)G - (long)A] * x2[(long)G - (long)A] +
	       freqt * x1[(long)T - (long)A] * x2[(long)T - (long)A];
      aa = prod12 - prod3;
      bb = prod3 - prod1*prod2;
      cc = prod1 * prod2;
      term[i][j] = zz * aa + z1 * bb + cc;
      slopeterm[i][j] = zzs * aa + z1s * bb;
      curveterm[i][j] = zzc * aa + z1c * bb;
    }
    sumterm = 0.0;
    for (j = 0; j < rcategs; j++)
      sumterm += probcat[j] * term[i][j];
    lterm = log(sumterm);
    for (j = 0; j < rcategs; j++) {
      term[i][j] = term[i][j] / sumterm;
      slopeterm[i][j] = slopeterm[i][j] / sumterm;
      curveterm[i][j] = curveterm[i][j] / sumterm; 
    }
    sum += aliasweight[i] * lterm;
  }
  for (i = 0; i < rcategs; i++) {
    thelike[i] = 1.0;
    theslope[i] = 0.0;
    thecurve[i] = 0.0;
  }
  for (i = 0; i < sites; i++) {
    sumc = 0.0;
    sumcs = 0.0;
    sumcc = 0.0;
    for (k = 0; k < rcategs; k++) {
      sumc += probcat[k] * thelike[k];
      sumcs += probcat[k] * theslope[k];
      sumcc += probcat[k] * thecurve[k];
    }
    sumc *= lambda;
    sumcs *= lambda;
    sumcc *= lambda;
    if ((ally[i] > 0) && (location[ally[i]-1] > 0)) {
      lai = location[ally[i] - 1];
      memcpy(clai, term[lai - 1], rcategs*sizeof(double));
      memcpy(cslai, slopeterm[lai - 1], rcategs*sizeof(double));
      memcpy(cclai, curveterm[lai - 1], rcategs*sizeof(double));
      if (weight[i] > 1) {
        for (j = 0; j < rcategs; j++) {
          if (clai[j] > 0.0)
            clai[j] = exp(weight[i]*log(clai[j]));
          else clai[j] = 0.0;
          if (cslai[j] > 0.0)
            cslai[j] = exp(weight[i]*log(cslai[j]));
          else cslai[j] = 0.0;
          if (cclai[j] > 0.0)
            cclai[j] = exp(weight[i]*log(cclai[j])); 
          else cclai[j] = 0.0;
      }
    }
      for (j = 0; j < rcategs; j++) {
        nulike[j] = ((1.0 - lambda) * thelike[j] + sumc) * clai[j];
        nuslope[j] = ((1.0 - lambda) * theslope[j] + sumcs) * clai[j]
                   + ((1.0 - lambda) * thelike[j] + sumc) * cslai[j];
        nucurve[j] = ((1.0 - lambda) * thecurve[j] + sumcc) * clai[j]
             + 2.0 * ((1.0 - lambda) * theslope[j] + sumcs) * cslai[j]
                   + ((1.0 - lambda) * thelike[j] + sumc) * cclai[j];
      }
    } else {
      for (j = 0; j < rcategs; j++) {
        nulike[j] = ((1.0 - lambda) * thelike[j] + sumc);
        nuslope[j] = ((1.0 - lambda) * theslope[j] + sumcs);
        nucurve[j] = ((1.0 - lambda) * thecurve[j] + sumcc);
      }
    }
    memcpy(thelike, nulike, rcategs*sizeof(double));
    memcpy(theslope, nuslope, rcategs*sizeof(double));
    memcpy(thecurve, nucurve, rcategs*sizeof(double));
  }
  sum2 = 0.0;
  slope2 = 0.0;
  curve2 = 0.0;
  for (i = 0; i < rcategs; i++) {
    sum2 += probcat[i] * thelike[i];
    slope2 += probcat[i] * theslope[i];
    curve2 += probcat[i] * thecurve[i];
  }
  sum += log(sum2);
  (*like) = sum;
  (*slope) = slope2 / sum2;
  (*curve) = (curve2 - slope2 * slope2 / sum2) / sum2;
} /* slopecurv */


void makenewv(node *p)
{
  /* Newton-Raphson algorithm improvement of a branch length */
  long it, ite;
  double y, yold=0, yorig, like, slope, curve, oldlike=0;
  boolean done, firsttime, better;
  node *q;

  q = p->back;
  y = p->v;
  yorig = y;
  done = false;
  firsttime = true;
  it = 1;
  ite = 0;
  while ((it < iterations) && (ite < 20) && (!done)) {
    slopecurv (p, y, &like, &slope, &curve);
    better = false;
    if (firsttime) {
      yold = y;
      oldlike = like;
      firsttime = false;
      better = true;
    } else {
      if (like > oldlike) {
        yold = y;
        oldlike = like;
        better = true;
        it++;
      }
    }
    if (better) {
      y = y + slope/fabs(curve);
      if (yold < epsilon)
        yold = epsilon;
    } else {
      if (fabs(y - yold) < epsilon)
        ite = 20;
      y = (y + 19*yold) / 20.0;
    }
    ite++;
    done = fabs(y-yold) < epsilon;
  }
  smoothed = (fabs(yold-yorig) < epsilon) && (yorig > 1000.0*epsilon);
  p->v = yold;
  q->v = yold;
  curtree.likelihood = oldlike;
}  /* makenewv */


void update(node *p)
{
  long num_sibs, i;
  node *sib_ptr;

  if (!p->tip && !p->initialized)
    nuview(p);
  if (!p->back->tip && !p->back->initialized)
    nuview(p->back);
  if ((!usertree) || (usertree && !lngths) || p->iter) {
    makenewv(p);

    if (!p->tip) {
      num_sibs = count_sibs (p);
      sib_ptr  = p;
      for (i=0; i < num_sibs; i++) {
	sib_ptr              = sib_ptr->next;
	sib_ptr->initialized = false;
      }
    }

    if (!p->back->tip) {
      num_sibs = count_sibs (p->back);
      sib_ptr  = p->back;
      for (i=0; i < num_sibs; i++) {
	sib_ptr              = sib_ptr->next;
        sib_ptr->initialized = false;
      }
    }
  }
}  /* update */


void smooth(node *p)
{
  long i, num_sibs;
  node *sib_ptr;
  
  smoothed = false;
  update (p);
  if (p->tip)
    return;

  num_sibs = count_sibs (p);
  sib_ptr  = p;
  
  for (i=0; i < num_sibs; i++) {
    sib_ptr = sib_ptr->next;
    
    if (polishing || (smoothit && !smoothed)) {
      smooth(sib_ptr->back);
      p->initialized = false;
      sib_ptr->initialized = false;
    }
  }
}  /* smooth */


void insert_(node *p, node *q, boolean dooinit)
{
  /* Insert q near p */
  long i, j, num_sibs;
  node *r, *sib_ptr;

  r = p->next->next;
  hookup(r, q->back);
  hookup(p->next, q);
  q->v = 0.5 * q->v;
  q->back->v = q->v;
  r->v = q->v;
  r->back->v = r->v;
  p->initialized = false;
  p->next->initialized = false;
  p->next->next->initialized = false;
  if (dooinit) {
    inittrav(p);
    inittrav(q);
    inittrav(q->back);
  }
  i = 1;
  while (i <= smoothings) {
    smooth (p);
    if (!smoothit) {
      if (!p->tip) {
        num_sibs = count_sibs (p);
        sib_ptr  = p;
        for (j=0; j < num_sibs; j++) {
          smooth (sib_ptr->next->back);
          sib_ptr = sib_ptr->next;
        }
      }
    }
    else 
      smooth(p->back);
    i++;
  }
}  /* insert_ */


void dnaml_re_move(node **p, node **q)
{
  /* remove p and record in q where it was */
  long i;

  /** assumes bifurcations */
  *q = (*p)->next->back;
  hookup(*q, (*p)->next->next->back);
  (*p)->next->back = NULL;
  (*p)->next->next->back = NULL;
  inittrav((*q));
  inittrav((*q)->back);
  i = 1;
  while (i <= smoothings) {
    smooth(*q);
    if (smoothit)
      smooth((*q)->back);
    i++;
  }
}  /* dnaml_re_move */


void buildnewtip(long m, tree *tr)
{
  node *p;

  p = tr->nodep[nextsp + spp - 3];
  hookup(tr->nodep[m - 1], p);
  p->v = initialv;
  p->back->v = initialv;
}  /* buildnewtip */


void buildsimpletree(tree *tr)
{
  hookup(tr->nodep[enterorder[0] - 1], tr->nodep[enterorder[1] - 1]);
  tr->nodep[enterorder[0] - 1]->v = 0.1;
  tr->nodep[enterorder[0] - 1]->back->v = 0.1;
  tr->nodep[enterorder[1] - 1]->v = 0.1;
  tr->nodep[enterorder[1] - 1]->back->v = 0.1;
  buildnewtip(enterorder[2], tr);
  insert_(tr->nodep[enterorder[2] - 1]->back,
          tr->nodep[enterorder[0] - 1], false);
}  /* buildsimpletree2 */


void addtraverse(node *p, node *q, boolean contin)
{
  /* try adding p at q, proceed recursively through tree */
  long i, num_sibs;
  double like, vsave = 0;
  node *qback = NULL, *sib_ptr;

  if (!smoothit) {
    vsave = q->v;
    qback = q->back;
  }
  insert_(p, q, false);
  like = evaluate(p, false);
  if (like > bestyet) {
    bestyet = like;
    if (smoothit)
      dnamlcopy(&curtree, &bestree, nonodes2, rcategs);
    else
      qwhere = q;
    succeeded = true;
  }
  if (smoothit)
    dnamlcopy(&priortree, &curtree, nonodes2, rcategs);
  else {
    hookup (q, qback);
    q->v = vsave;
    q->back->v = vsave;
    curtree.likelihood = bestyet;
  }
  if (!q->tip && contin) {
    num_sibs = count_sibs (q);
    if (q == curtree.start)
      num_sibs++;
    sib_ptr  = q;
    for (i=0; i < num_sibs; i++) {
      addtraverse(p, sib_ptr->next->back, contin);
      sib_ptr = sib_ptr->next;
    }
  }

}  /* addtraverse */


void rearrange(node *p, node *pp)
{
  /* rearranges the tree, globally or locally moving pp around near p */
  long i, num_sibs;
  double v3 = 0, v4 = 0, v5 = 0;
  node *q, *r, *sib_ptr;

  if (!p->tip && !p->back->tip) {
    curtree.likelihood = bestyet;
    if (p->back->next != pp)
      r = p->back->next;
    else
      r = p->back->next->next;
    /* assumes bifurcations? */
    if (!smoothit) {
      v3 = r->v;
      v4 = r->next->v;
      v5 = r->next->next->v;
    }
    else
      dnamlcopy(&curtree, &bestree, nonodes2, rcategs);
    dnaml_re_move(&r, &q);
    if (smoothit)
      dnamlcopy(&curtree, &priortree, nonodes2, rcategs);
    else
      qwhere = q;
    num_sibs = count_sibs (p);
    sib_ptr  = p;
    for (i=0; i < num_sibs; i++) {
      sib_ptr = sib_ptr->next;
      addtraverse(r, sib_ptr->back, global && nextsp == spp);
    }
    if (global && nextsp == spp && !succeeded) {
      p = p->back;
      if (!p->tip) {
        num_sibs = count_sibs (p);
        sib_ptr  = p;
        for (i=0; i < num_sibs; i++) {
          sib_ptr = sib_ptr->next;
          addtraverse(r, sib_ptr->back, global && nextsp == spp);
        }
      }
      p = p->back;
    }
    if (smoothit)
      dnamlcopy(&bestree, &curtree, nonodes2, rcategs);
    else {
      insert_(r, qwhere, true);
      if (qwhere == q) {
        r->v = v3;
        r->back->v = v3;
        r->next->v = v4;
        r->next->back->v = v4;
        r->next->next->v = v5;
        r->next->next->back->v = v5;
        curtree.likelihood = bestyet;
      }
      else {
        smoothit = true;
        for (i = 1; i<=smoothings; i++) {
          smooth (r);
          smooth (r->back);
        }
        smoothit = false;
        dnamlcopy(&curtree, &bestree, nonodes2, rcategs);
      }
    }
    if (global && nextsp == spp && progress) {
      putchar('.');
      fflush(stdout);
    }
  }
  if (!p->tip) {
    num_sibs = count_sibs (p);
    if (p == curtree.start)
      num_sibs++;
    sib_ptr  = p;
    for (i=0; i < num_sibs; i++) {
      sib_ptr = sib_ptr->next;
      rearrange(sib_ptr->back, p);
    }
  }
}  /* rearrange */


void initdnamlnode(node **p, node **grbg, node *q, long len, long nodei,
			long *ntips, long *parens, initops whichinit,
			pointarray treenode, pointarray nodep, Char *str, Char *ch,
			FILE *intree)
{
  /* initializes a node */
  boolean minusread;
  double valyew, divisor;

  switch (whichinit) {
  case bottom:
    gnu(grbg, p);
    (*p)->index = nodei;
    (*p)->tip = false;
    malloc_pheno((*p), spp, endsite, rcategs);
    nodep[(*p)->index - 1] = (*p);
    break;
  case nonbottom:
    gnu(grbg, p);
    malloc_pheno(*p, spp, endsite, rcategs);
    (*p)->index = nodei;
    break;
  case tip:
    match_names_to_data (str, nodep, p, spp);
    break;
  case iter:
    (*p)->initialized = false;
    (*p)->v = initialv;
    (*p)->iter = true;
    if ((*p)->back != NULL)
      (*p)->back->iter = true;
    break;
  case length:
    processlength(&valyew, &divisor, ch, &minusread, intree, parens);
    (*p)->v = valyew / divisor / fracchange;
    (*p)->iter = false;
    if ((*p)->back != NULL) {
      (*p)->back->v = (*p)->v;
      (*p)->back->iter = false;
    }
    break;
  case hsnolength:
    haslengths = false;
    break;
  default:	/* cases hslength, treewt, unittrwt */
    break;	/* should never occur				*/
  }
} /* initdnamlnode */


void dnaml_coordinates(node *p, double lengthsum, long *tipy,
			double *tipmax)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  double xx;

  if (p->tip) {
    p->xcoord = (long)(over * lengthsum + 0.5);
    p->ycoord = (*tipy);
    p->ymin = (*tipy);
    p->ymax = (*tipy);
    (*tipy) += down;
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;
  do {
    xx = fracchange * q->v;
    if (xx > 100.0)
      xx = 100.0;
    dnaml_coordinates(q->back, lengthsum + xx, tipy,tipmax);
    q = q->next;
  } while ((p == curtree.start || p != q) &&
	   (p != curtree.start || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (long)(over * lengthsum + 0.5);
  if (p == curtree.start)
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* dnaml_coordinates */


void dnaml_printree()
{
  /* prints out diagram of the tree2 */
  long tipy;
  double scale, tipmax;
  long i;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  dnaml_coordinates(curtree.start, 0.0, &tipy, &tipmax);
  scale = 1.0 / (long)(tipmax + 1.000);
  for (i = 1; i <= (tipy - down); i++)
    drawline2(i, scale, curtree);
  putc('\n', outfile);
}  /* dnaml_printree */


void sigma(node *p, double *sumlr, double *s1, double *s2)
{
  /* compute standard deviation */
  double tt, aa, like, slope, curv;

  slopecurv (p, p->v, &like, &slope, &curv);
  tt = p->v;
  p->v = epsilon;
  p->back->v = epsilon;
  aa = evaluate(p, false);
  p->v = tt;
  p->back->v = tt;
  (*sumlr) = evaluate(p, false) - aa;
  if (curv < -epsilon) {
    (*s1) = p->v + (-slope - sqrt(slope * slope -  3.841 * curv)) / curv;
    (*s2) = p->v + (-slope + sqrt(slope * slope -  3.841 * curv)) / curv;
  }
  else {
    (*s1) = -1.0;
    (*s2) = -1.0;
  }
}  /* sigma */


void describe(node *p)
{
  /* print out information for one branch */
  long i, num_sibs;
  node *q, *sib_ptr;
  double sumlr, sigma1, sigma2;

  if (!p->tip && !p->initialized)
    nuview(p);
  if (!p->back->tip && !p->back->initialized)
    nuview(p->back);
  q = p->back;
  if (q->tip) {
    fprintf(outfile, " ");
    for (i = 0; i < nmlngth; i++)
      putc(nayme[q->index-1][i], outfile);
    fprintf(outfile, "    ");
  } else
    fprintf(outfile, "  %4ld          ", q->index - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index-1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "%15.5f", q->v * fracchange);
  if (!usertree || (usertree && !lngths) || p->iter) {
    sigma(q, &sumlr, &sigma1, &sigma2);
    if (sigma1 <= sigma2)
      fprintf(outfile, "     (     zero,    infinity)");
    else {
      fprintf(outfile, "     (");
      if (sigma2 <= 0.0)
        fprintf(outfile, "     zero");
      else
        fprintf(outfile, "%9.5f", sigma2 * fracchange);
      fprintf(outfile, ",%12.5f", sigma1 * fracchange);
      putc(')', outfile);
      }
    if (sumlr > 1.9205)
      fprintf(outfile, " *");
    if (sumlr > 2.995)
      putc('*', outfile);
    }
  putc('\n', outfile);
  if (!p->tip) {
    num_sibs = count_sibs (p);
    sib_ptr  = p;
    for (i=0; i < num_sibs; i++) {
      sib_ptr = sib_ptr->next;
      describe(sib_ptr->back);
    }
  }
}  /* describe */


void reconstr(node *p, long n)
{
  /* reconstruct and print out base at site n+1 at node p */
  long i, j, m, first, second;
  double f, sum, xx[4];

  if ((ally[n] == 0) || (location[ally[n]-1] == 0))
    putc('.', outfile);
  else {
    j = location[ally[n]-1] - 1;
    for (i = 0; i < 4; i++) {
      f = p->x[j][mx-1][i];
      if (!p->tip)
        f = f * p->next->x[j][mx-1][i] * p->next->next->x[j][mx-1][i];
      f = sqrt(f);
      xx[i] = f;
    }
    xx[0] *= freqa;
    xx[1] *= freqc;
    xx[2] *= freqg;
    xx[3] *= freqt;
    sum = xx[0]+xx[1]+xx[2]+xx[3];
    for (i = 0; i < 4; i++)
      xx[i] /= sum;
    first = 0;
    for (i = 1; i < 4; i++)
      if (xx [i] > xx[first])
        first = i;
    if (first == 0)
      second = 1;
    else
      second = 0;
    for (i = 0; i < 4; i++)
      if ((i != first) && (xx[i] > xx[second]))
        second = i;
    m = 1 << first;
    if (xx[first] < 0.4999995)
      m = m + (1 << second);
    if (xx[first] > 0.95)
      putc(toupper(basechar[m - 1]), outfile);
    else
      putc(basechar[m - 1], outfile);
    if (rctgry && rcategs > 1)
      mx = mp[n - 1][mx - 1];    
    else
      mx = 1;
  }
} /* reconstr */


void rectrav(node *p, long m, long n)
{
  /* print out segment of reconstructed sequence for one branch */
  long i;

  putc(' ', outfile);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index-1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "  ");
  mx = mx0;
  for (i = m; i <= n; i++) {
    if ((i % 10 == 0) && (i != m))
      putc(' ', outfile);
    reconstr(p, i);
  }
  putc('\n', outfile);
  if (!p->tip) {
    rectrav(p->next->back, m, n);
    rectrav(p->next->next->back, m, n);
  }
}  /* rectrav */


void summarize()
{
  /* print out branch length information and node numbers */
  long i, j, num_sibs;
  double mode, sum;
  double like[maxcategs],nulike[maxcategs];
  double **marginal;
  node   *sib_ptr;

  if (!treeprint)
    return;
  fprintf(outfile, "\nremember: ");
  if (outgropt)
    fprintf(outfile, "(although rooted by outgroup) ");
  fprintf(outfile, "this is an unrooted tree!\n\n");
  fprintf(outfile, "Ln Likelihood = %11.5f\n", curtree.likelihood);
  fprintf(outfile, "\n Between        And            Length");
  if (!(usertree && lngths && haslengths))
    fprintf(outfile, "      Approx. Confidence Limits");
  fprintf(outfile, "\n");
  fprintf(outfile, " -------        ---            ------");
  if (!(usertree && lngths && haslengths))
    fprintf(outfile, "      ------- ---------- ------");
  fprintf(outfile, "\n\n");
  for (i = spp; i < nonodes2; i++) {
    /* So this works with arbitrary multifurcations */
    if (curtree.nodep[i]) {
      num_sibs = count_sibs (curtree.nodep[i]);
      sib_ptr  = curtree.nodep[i];
      for (j = 0; j < num_sibs; j++) {
        sib_ptr->initialized = false;
        sib_ptr              = sib_ptr->next;
      }
    }
  }

  describe(curtree.start->back);

  /* So this works with arbitrary multifurcations */
  num_sibs = count_sibs (curtree.start);
  sib_ptr  = curtree.start;
  for (i=0; i < num_sibs; i++) {
    sib_ptr = sib_ptr->next;
    describe(sib_ptr->back);
  }

  fprintf(outfile, "\n");
  if (!(usertree && lngths && haslengths)) {
    fprintf(outfile, "     *  = significantly positive, P < 0.05\n");
    fprintf(outfile, "     ** = significantly positive, P < 0.01\n\n");
  }
  dummy = evaluate(curtree.start, false);
  if (rctgry && rcategs > 1) {
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++) {
	nulike[j] = (1.0 - lambda + lambda * probcat[j]) * like[j];
	mp[i][j] = j + 1;
	for (k = 1; k <= rcategs; k++) {
	  if (k != j + 1) {
	    if (lambda * probcat[k - 1] * like[k - 1] > nulike[j]) {
	      nulike[j] = lambda * probcat[k - 1] * like[k - 1];
	      mp[i][j] = k;
	    }
	  }
	}
        if ((ally[i] > 0) && (location[ally[i]-1] > 0))
	  nulike[j] *= contribution[location[ally[i] - 1] - 1][j];
	sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
	nulike[j] /= sum;
      memcpy(like, nulike, rcategs * sizeof(double));
    }
    mode = 0.0;
    mx = 1;
    for (i = 1; i <= rcategs; i++) {
      if (probcat[i - 1] * like[i - 1] > mode) {
	mx = i;
	mode = probcat[i - 1] * like[i - 1];
      }
    }
    mx0 = mx;
    fprintf(outfile,
 "Combination of categories that contributes the most to the likelihood:\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', outfile);
    for (i = 1; i <= sites; i++) {
      fprintf(outfile, "%ld", mx);
      if (i % 10 == 0)
	putc(' ', outfile);
      if (i % 60 == 0 && i != sites) {
	putc('\n', outfile);
	for (j = 1; j <= nmlngth + 3; j++)
	  putc(' ', outfile);
      }
      mx = mp[i - 1][mx - 1];
    }
    fprintf(outfile, "\n\n");
    marginal = (double **) calloc (1,sites*sizeof(double *));
    for (i = 0; i < sites; i++)
      marginal[i] = (double *) calloc (1,rcategs*sizeof(double));
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++) {
	nulike[j] = (1.0 - lambda + lambda * probcat[j]) * like[j];
	for (k = 1; k <= rcategs; k++) {
	  if (k != j + 1)
	      nulike[j] += lambda * probcat[k - 1] * like[k - 1];
	}
        if ((ally[i] > 0) && (location[ally[i]-1] > 0))
	  nulike[j] *= contribution[location[ally[i] - 1] - 1][j];
	sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++) {
	nulike[j] /= sum;
        marginal[i][j] = nulike[j];
      }
      memcpy(like, nulike, rcategs * sizeof(double));
    }
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = 0; i < sites; i++) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++) {
	nulike[j] = (1.0 - lambda + lambda * probcat[j]) * like[j];
	for (k = 1; k <= rcategs; k++) {
	  if (k != j + 1)
	      nulike[j] += lambda * probcat[k - 1] * like[k - 1];
	}
        marginal[i][j] *= like[j] * probcat[j];
	sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
	nulike[j] /= sum;
      memcpy(like, nulike, rcategs * sizeof(double));
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
        sum += marginal[i][j];
      for (j = 0; j < rcategs; j++)
        marginal[i][j] /= sum;
    }
    fprintf(outfile, "Most probable category at each site if > 0.95 probability (\".\" otherwise)\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', outfile);
    for (i = 0; i < sites; i++) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
        if (marginal[i][j] > sum) {
          sum = marginal[i][j];
          mx0 = j;
        }
        if (sum >= 0.95)
        fprintf(outfile, "%ld", mx0+1);
      else
        putc('.', outfile);
      if ((i+1) % 60 == 0) {
        if (i != 0) {
	  putc('\n', outfile);
	  for (j = 1; j <= nmlngth + 3; j++)
	    putc(' ', outfile);
        }
      }
      else if ((i+1) % 10 == 0)
	putc(' ', outfile);
    }
    putc('\n', outfile);
    for (i = 0; i < sites; i++)
      free(marginal[i]);
    free(marginal);
  }
  putc('\n', outfile);
  if (hypstate) {
    fprintf(outfile, "Probable sequences at interior nodes:\n\n");
    fprintf(outfile, "  node      ");
    for (i = 0; (i < 14) || (i < ((sites-1) % 2)-19) ; i++)
      putc(' ', outfile);
    fprintf(outfile, "Reconstructed sequence (caps if > 0.95)\n\n");
    if (!rctgry || (rcategs == 1))
      mx0 = 1;
    for (i = 0; i < sites; i += 60) {
      k = i + 59;
      if (k >= sites)
        k = sites - 1;
      rectrav(curtree.start, i, k);
      rectrav(curtree.start->back, i, k);
      putc('\n', outfile);
    }
  }
}  /* summarize */


void dnaml_treeout(node *p)
{
  /* write out file with representation of final tree2 */
  /* Only works for bifurcations! */
  long i, n, w;
  Char c;
  double x;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[p->index-1][i - 1] != ' ')
	n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index-1][i];
      if (c == ' ')
	c = '_';
      putc(c, outtree);
    }
    col += n;
  } else {
    putc('(', outtree);
    col++;
    dnaml_treeout(p->next->back);
    putc(',', outtree);
    col++;
    if (col > 45) {
      putc('\n', outtree);
      col = 0;
    }
    dnaml_treeout(p->next->next->back);
    if (p == curtree.start) {
      putc(',', outtree);
      col++;
      if (col > 45) {
	putc('\n', outtree);
	col = 0;
      }
      dnaml_treeout(p->back);
    }
    putc(')', outtree);
    col++;
  }
  x = p->v * fracchange;
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p == curtree.start)
    fprintf(outtree, ";\n");
  else {
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    col += w + 8;
  }
}  /* dnaml_treeout */


void inittravtree(node *p)
{
  /* traverse tree to set initialized and v to initial values */

  p->initialized = false;
  p->back->initialized = false;
  if (!p->tip) {
    inittravtree(p->next->back);
    inittravtree(p->next->next->back);
  }
} /* inittravtree */


void treevaluate()
{
  /* evaluate a user tree */
  long i;

  inittravtree(curtree.start);
  polishing = true;
  smoothit = true;
  for (i = 1; i <= smoothings * 4; i++)
    smooth (curtree.start);
  dummy = evaluate(curtree.start, true);
}  /* treevaluate */


void maketree()
{
  long i, j;
  boolean dummy_first, goteof;
  pointarray dummy_treenode;
  long nextnode;
  node *root, *q, *r;

  inittable();

  if (usertree) {
    openfile(&intree,INTREE,"input tree file", "r",progname,intreename);
    inittable_for_usertree (intree, progname, intreename);
    numtrees = countsemic(&intree);
    l0gl = (double *) calloc (1,numtrees * sizeof(double));
    l0gf = (double **) calloc (1,numtrees * sizeof(double *));
    for (i=0; i < numtrees; ++i)
      l0gf[i] = (double *) calloc (1,endsite * sizeof(double));
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
	putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    which = 1;

    /* This taken out of tree read, used to be [spp-1], but referring
       to [0] produces output identical to what the pre-modified dnaml
       produced. */

    while (which <= numtrees) {
      
      /* These initializations required each time through the loop
         since multiple trees require re-initialization */
      haslengths = true;
      nextnode         = 0;
      dummy_first      = true;
      goteof           = false;

      treeread(intree, &root, dummy_treenode, &goteof, &dummy_first,
               curtree.nodep, INTREE, "Dnaml", &nextnode,
               &haslengths, &grbg, initdnamlnode);
      q = root;
      r = root;
      while (!(q->next == root))
        q = q->next;
      q->next = root->next;
      root = q;
      chuck(&grbg, r);
      curtree.nodep[spp] = q;
      if (goteof && (which <= numtrees)) {
	/* if we hit the end of the file prematurely */
	printf ("\n");
	printf ("ERROR: trees missing at end of file.\n");
	printf ("\tExpected number of trees:\t\t%ld\n", numtrees);
	printf ("\tNumber of trees actually in file:\t%ld.\n\n", which - 1);
	exxit(-1);
      }

      curtree.start = curtree.nodep[0]->back;
      
      treevaluate();
      if (reconsider) {
        bestyet = nextsp*sites*(freqa*log(freqa)+freqg*log(freqg) +
				freqc*log(freqc)+freqt*log(freqt));
        succeeded = true;
        while (succeeded) {
          succeeded = false;
          rearrange(curtree.start, curtree.start->back);
        }
        treevaluate();
      }
      dnaml_printree();
      summarize();
      if (trout) {
	col = 0;
	dnaml_treeout(curtree.start);
      }
      which++;
    }
    FClose(intree);
    putc('\n', outfile);
    if (!auto_ && numtrees > 1 && weightsum > 1 )
      standev2(numtrees, weightsum, maxwhich, 0, endsite-1, maxlogl,
               l0gl, l0gf, aliasweight);
  } else {
    /* If there's no input user tree, */
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    if (progress) {
      printf("\nAdding species:\n");
      writename(0, 3, enterorder);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    nextsp = 3;
    polishing = false;
    buildsimpletree(&curtree);
    curtree.start = curtree.nodep[enterorder[0] - 1]->back;
    smoothit = improve;
    nextsp = 4;
    while (nextsp <= spp) {
      buildnewtip(enterorder[nextsp - 1], &curtree);
      bestyet = nextsp*sites*(freqa*log(freqa)+freqg*log(freqg)
                                +freqc*log(freqc)+freqt*log(freqt));
      if (smoothit)
        dnamlcopy(&curtree, &priortree, nonodes2, rcategs);
      addtraverse(curtree.nodep[enterorder[nextsp - 1] - 1]->back,
		  curtree.start, true);
      if (smoothit)
        dnamlcopy(&bestree, &curtree, nonodes2, rcategs);
      else {
        insert_(curtree.nodep[enterorder[nextsp - 1] - 1]->back, qwhere, true);
        smoothit = true;
        for (i = 1; i<=smoothings; i++) {
          smooth (curtree.start);
          smooth (curtree.start->back);
        }
        smoothit = false;
        dnamlcopy(&curtree, &bestree, nonodes2, rcategs);
        bestyet = curtree.likelihood;
      }
      if (progress) {
        writename(nextsp - 1, 1, enterorder);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
      if (global && nextsp == spp && progress) {
	printf("Doing global rearrangements\n");
        printf("  !");
        for (j = 1; j <= (spp - 3); j++)
          putchar('-');
        printf("!\n");
      }
      succeeded = true;
      while (succeeded) {
	succeeded = false;
        if (global && nextsp == spp && progress) {
          printf("   ");
          fflush(stdout);
        }
        rearrange(curtree.start, curtree.start->back);
        if (global && nextsp == spp && progress)
          putchar('\n');
      }
      for (i = spp; i < nextsp + spp - 2; i++) {
        curtree.nodep[i]->initialized = false;
        curtree.nodep[i]->next->initialized = false;
        curtree.nodep[i]->next->next->initialized = false;
      }
      if (!smoothit) {
        smoothit = true;
        for (i = 1; i<=smoothings; i++) {
          smooth (curtree.start);
          smooth (curtree.start->back);
        }
        smoothit = false;
        dnamlcopy(&curtree, &bestree, nonodes2, rcategs);
        bestyet = curtree.likelihood;
      }
      nextsp++;
    }
    if (global && progress) {
      putchar('\n');
      fflush(stdout);
    }
    if (njumble > 1) {
      if (jumb == 1)
        dnamlcopy(&bestree, &bestree2, nonodes2, rcategs);
      else
        if (bestree2.likelihood < bestree.likelihood)
          dnamlcopy(&bestree, &bestree2, nonodes2, rcategs);
    }
    if (jumb == njumble) {
      if (njumble > 1)
        dnamlcopy(&bestree2, &curtree, nonodes2, rcategs);
      curtree.start = curtree.nodep[outgrno - 1]->back;
      for (i = 0; i < nonodes2; i++) {
        if (i < spp)
          curtree.nodep[i]->initialized = false;
        else {
          curtree.nodep[i]->initialized = false;
          curtree.nodep[i]->next->initialized = false;
          curtree.nodep[i]->next->next->initialized = false;
        } 
      }
      treevaluate();
      dnaml_printree();
      summarize();
      if (trout) {
	col = 0;
	dnaml_treeout(curtree.start);
      }
    }
  }
  if (usertree) {
    free(l0gl);
    for (i=0; i < numtrees; i++)
      free(l0gf[i]);
    free(l0gf);
  }
  for (i = 0; i < rcategs; i++)
    for (j = 0; j < categs; j++)
      free(tbl[i][j]);
  for (i = 0; i < rcategs; i++)
    free(tbl[i]);
  free(tbl);
  if (jumb < njumble)
    return;
  free(contribution);
  free(mp);
  for (i=0; i < endsite; i++)
     free(term[i]);
  free(term);
  for (i=0; i < endsite; i++)
     free(slopeterm[i]);
  free(slopeterm);
  for (i=0; i < endsite; i++)
     free(curveterm[i]);
  free(curveterm);
  free_all_x_in_array (nonodes2, curtree.nodep);
  if (!usertree || reconsider) {
    free_all_x_in_array (nonodes2, bestree.nodep);
    free_all_x_in_array (nonodes2, priortree.nodep);
    if (njumble > 1)
      free_all_x_in_array (nonodes2, bestree2.nodep);
  }
  if (progress) {
    printf("\n\nOutput written to output file\n\n");
    if (trout)
      printf("Tree also written onto file\n");
    putchar('\n');
  }
}  /* maketree */


void clean_up()
{
  /* Free and/or close stuff */
  long i;

    free (rrate);
    free (probcat);
    free (rate);
    /* Seems to require freeing every time... */
    for (i = 0; i < spp; i++) {
      free (y[i]);
    }
    free (y);
    free (nayme);
    free (enterorder);
    free (category);
    free (weight);
    free (alias);
    free (ally);
    free (location);
    free (aliasweight);

#if 0          /* ???? debug ???? */
  freetree2(curtree.nodep, nonodes2);

  if (! (usertree && !reconsider)) {
    freetree2(bestree.nodep, nonodes2);
    freetree2(priortree.nodep, nonodes2);
  }
  
  if (! (njumble <= 1))
    freetree2(bestree2.nodep, nonodes2);
#endif
  FClose(infile);
  FClose(outfile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
}   /* clean_up */


int main(int argc, Char *argv[])
{  /* DNA Maximum Likelihood */
#ifdef MAC
  argc = 1;		/* macsetup("DnaML","");	*/
  argv[0] = "DnaML";
#endif
#ifdef WIN32
  phySetConsoleAttributes();
  phyClearScreen();
#endif
  strcpy(progname,argv[0]);
  openfile(&infile,INFILE,"input file","r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file","w",argv[0],outfilename);
  mulsets = false;
  datasets = 1;
  firstset = true;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  grbg = NULL;
  doinit();
  ttratio0 = ttratio;
  if (ctgry)
    openfile(&catfile,CATFILE,"categories file","r",argv[0],catfilename);
  if (weights || justwts)
    openfile(&weightfile,WEIGHTFILE,"weights file","r",argv[0],weightfilename);
  if (trout)
    openfile(&outtree,OUTTREE,"output tree file","w",argv[0],outtreename);
  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n", ith);
      printf("\nData set # %ld:\n", ith);
    }
    ttratio = ttratio0;
    getinput();
    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
  }

  clean_up();
  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* DNA Maximum Likelihood */
 
