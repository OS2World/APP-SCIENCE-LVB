#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1986-2000 by the University of Washington
  and by Joseph Felsenstein.  Written by Joseph Felsenstein.  Permission is
  granted to copy and use this program provided no fee is charged for it
  and provided that this copyright notice is not removed. */

#define epsilon         0.0001   /* used in makenewv, getthree, update */
#define over            60

typedef struct valrec {
  double rat, ratxi, ratxv, orig_zz, z1, y1, z1zz, z1yy, xiz1,
         xiy1xv;
  double *ww, *zz, *wwzz, *vvzz; 
} valrec;

typedef long vall[maxcategs];

typedef double contribarr[maxcategs];

extern FILE *infile, *outfile, *intree, *outtree;
extern long spp, nonodes, endsite;
extern boolean interleaved, printdata, treeprint;
extern sequence y;
extern steptr weight, alias, location, ally;
extern naym *nayme;
valrec ***tbl;

extern void initcategs() ;
extern void justweights();
extern void inputcategs();
extern void printcategs();
extern long countcomma();
extern void malloc_pheno();
extern void hookup();
extern void match_names_to_data();
extern long countsemic();
extern long count_sibs();
extern void freex();
extern void inittrav();
 
#ifndef OLDC
/* function prototypes */
boolean letter(Char);
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   inputoptions(void);
void   makeweights(void);
void   getinput(void);
void   inittable_for_usertree (FILE *, char *, char *);
void   inittable(void);
void   exmake(double, long);
void   alloc_nvd(long, nuview_data *);
void   free_nvd(long, nuview_data *);
void   nuview(node *);
double evaluate(node *);
void   oldnuview(node *);
void   getthree(node *);
void   makenewv(node *);
void   update(node *);
void   smooth(node *);
void   restoradd(node *, node *, node *, double);

void   dnamlk_add(node *, node *, node *, boolean );
void   dnamlk_re_move(node **, node **, boolean);
void   tryadd(node *, node **, node **);
void   addpreorder(node *, node *, node *, boolean, boolean);
void   tryrearr(node *, boolean *);
void   repreorder(node *, boolean *);
void   rearrange(node **);
void   initdnamlnode(node **, node **, node *, long, long, long *, long *,
		initops, pointarray, pointarray, Char *, Char *, FILE *);
void   tymetrav(node *, double *);
void   dnamlk_coordinates(node *, long *);

void   dnamlk_drawline(long, double);
void   dnamlk_printree(void);
void   describe(node *);
void   oldsummarize(void);
void   summarize(void);
void   dnamlk_treeout(node *);
void   nodeinit(node *);
void   oldnodeinit(node *);
void   initrav(node *);
void   travinit(node *);

void   travsp(node *);
void   treevaluate(void);
void   maketree(void);
/* function prototypes */
#endif


Char infilename[100], outfilename[100], intreename[100], outtreename[100],
     catfilename[100], weightfilename[100];
double *rrate;
long sites, weightsum, categs, datasets, ith, njumble, jumb, numtrees;
/*  sites = number of sites in actual sequences
  numtrees = number of user-defined trees */
long inseed;
extern boolean ibmpc, ansi;
boolean freqsfrom, global, jumble, lengths, trout, usertree, weights, rctgry,
               ctgry, ttr, auto_, progress, mulsets, firstset, reconsider,
               smoothit, polishing, justwts;
tree curtree, bestree, bestree2, priortree;
node *qwhere, *grbg;
double xi, xv, ttratio, ttratio0, freqa, freqc, freqg, freqt, freqr,
              freqy, freqar, freqcy, freqgr, freqty, fracchange, sumrates,
              lambda, lambda1;
long *enterorder;
steptr aliasweight;
double *rate;
double **term, **slopeterm, **curveterm;
longer seed;
double *probcat;
contribarr *contribution;
Char progname[20];
long rcategs, nonodes2;
vall *mp;


/* Local variables for maketree, propagated globally for C version: */
long    k, maxwhich, col;
double  like, bestyet, tdelta, lnlike, slope, curv, maxlogl;
boolean lastsp, smoothed, succeeded;
double  *l0gl;
double  x[3], lnl[3];
double  expon1i[maxcategs], expon1v[maxcategs],
        expon2i[maxcategs], expon2v[maxcategs];
node   *there;
double  **l0gf;
Char ch, ch2;


boolean letter(Char ch)
{
  /* tests whether ch is a letter
   -- works for both ASCII and EBCDIC codes */
  return ((ch >= 'A' && ch <= 'I') || (ch >= 'J' && ch <= 'R') ||
          (ch >= 'S' && ch <= 'Z') || (ch >= 'a' && ch <= 'i') ||
          (ch >= 'j' && ch <= 'r') || (ch >= 's' && ch <= 'z'));
}  /* letter */


void getoptions()
{
  /* interactively set options */
  long inseed0;
  Char ch;
  boolean done;
  boolean didchangecat, didchangercat;
  double probsum;

  fprintf(outfile, "\nNucleic acid sequence\n");
  fprintf(outfile, "   Maximum Likelihood method with molecular ");
  fprintf(outfile, "clock, version %s\n\n", VERSION);
  putchar('\n');
  auto_ = false;
  ctgry = false;
  didchangecat = false;
  rctgry = false;
  didchangercat = false;
  categs = 1;
  rcategs = 1;
  freqsfrom = true;
  global = false;
  jumble = false;
  njumble = 1;
  lambda = 1.0;
  lambda1 = 0.0;
  lengths = false;
  trout = true;
  ttratio = 2.0;
  ttr = false;
  usertree = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  interleaved = true;
  do {
    if (ansi)
      printf("\033[2J\033[H");
    else
      putchar('\n');
    printf("\nNucleic acid sequence\n");
    printf("   Maximum Likelihood method with molecular clock, version %s\n\n",
           VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?");
    if (usertree)
     printf("  No, use user trees in input file\n");
    else
      printf("  Yes\n");
    if (usertree) {
      printf("  L           Use lengths from user tree?");
      if (lengths)
        printf("  Yes\n");
      else
        printf("  No\n");
    }
    printf("  T        Transition/transversion ratio:");
    if (!ttr)
      printf("  2.0\n");
    else
      printf("  %8.4f\n", ttratio);
    printf("  F       Use empirical base frequencies?");
    if (freqsfrom)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  C   One category of substitution rates?");
    if (!ctgry)
      printf("  Yes\n");
    else
      printf("  %ld categories\n", categs);
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
    if (!usertree) {
      printf("  G                Global rearrangements?");
      if (global)
        printf("  Yes\n");
      else
        printf("  No\n");
      printf("  J   Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed = %8ld, %3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", datasets,
               (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?");
    if (interleaved)
      printf("  Yes\n");
    else
      printf("  No, sequential\n");
    printf("  0   Terminal type (IBM PC, ANSI, none)?");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi)
      printf("  ANSI\n");
    if (!(ibmpc || ansi))
      printf("  (none)\n");
    printf("  1    Print out the data at start of run");
    if (printdata)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  2  Print indications of progress of run");
    if (progress)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  3                        Print out tree");
    if (treeprint)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  4       Write out trees onto tree file?");
    if (trout)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      uppercase(&ch);
      if (strchr("JUCRAFGLTMI01234", ch) != NULL){
          switch (ch) {

        case 'C':
          ctgry = !ctgry;
          if (ctgry) {
            printf("\nSitewise user-assigned categories:\n\n");
            initcatn(&categs);
            if (rate){
              free(rate);
            }
            rate    = (double *) calloc (1, categs * sizeof(double));
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
            probcat = (double *) calloc (1, rcategs * sizeof(double));
            rrate   = (double *) calloc (1, rcategs * sizeof(double));
            didchangercat = true;
            initcategs(rcategs, rrate);
            initprobcat(rcategs, &probsum, probcat);
          }
          break;

        case 'A':
          auto_ = !auto_;
          if (auto_) {
            initlambda(&lambda);
            lambda1 = 1.0 - lambda;
          }
          break;

        case 'F':
          freqsfrom = !freqsfrom;
          if (!freqsfrom) {
              initfreqs(&freqa, &freqc, &freqg, &freqt);
          }
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
          lengths = !lengths;
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

        case 'M':
        mulsets = !mulsets;
        if (mulsets) {
          printf("Multiple data sets or multiple weights?");
          do {
            printf(" (type D or W)\n");
#ifdef WIN32
            phyFillScreenColor();
#endif
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
        }
      } else
        printf("Not a possible option!\n");
    }
  } while (!done);
  if (!didchangercat){
    rrate      = calloc (1, rcategs*sizeof(double));
    probcat    = calloc (1, rcategs*sizeof(double));
    rrate[0]   = 1.0;
    probcat[0] = 1.0;
  }
  if (!didchangecat){
    rate       = calloc (1, categs*sizeof(double));
    rate[0]    = 1.0;
  }
}  /* getoptions */


void allocrest()
{
  long i;

  y     = (Char **)Malloc(spp*sizeof(Char *));
  nayme  = (naym *)Malloc(spp*sizeof(naym));
  for (i = 0; i < spp; i++)
    y[i] = (char *)Malloc(sites * sizeof(char));
  enterorder  = (long *)Malloc(spp*sizeof(long));
  weight      = (long *)Malloc(sites*sizeof(long));
  category    = (long *)Malloc(sites*sizeof(long));
  alias       = (long *)Malloc(sites*sizeof(long));
  aliasweight = (long *)Malloc(sites*sizeof(long));
  ally        = (long *)Malloc(sites*sizeof(long));
  location    = (long *)Malloc(sites*sizeof(long));
}  /* allocrest */


void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &sites, &nonodes, 1);
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, sites);
  alloctree(&curtree.nodep, nonodes, usertree);
  allocrest();
  if (usertree)
    return;
  alloctree(&bestree.nodep, nonodes, 0);
  if (njumble <= 1)
    return;
  alloctree(&bestree2.nodep, nonodes, 0);
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
    inputcategs(0, sites, category, categs, "DnaMLK");
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
  sitesort2(sites, aliasweight);
  sitecombine2(sites, aliasweight);
  sitescrunch2(sites, 1, 2, aliasweight);
  for (i = 1; i <= sites; i++) {
    if (aliasweight[i - 1] > 0)
      endsite = i;
  }
  for (i = 1; i <= endsite; i++) {
    ally[alias[i - 1] - 1] = alias[i - 1];
    location[alias[i - 1] - 1] = i;
  }
  sumrates = 0.0;
  for (i = 0; i < rcategs; i++)
    sumrates += probcat[i] * rate[i];
  for (i = 0; i < rcategs; i++)
    rate[i] /= sumrates;
  contribution = (contribarr *) calloc (1, endsite*sizeof(contribarr));
}  /* makeweights */


void getinput()
{
  long grcategs;
  
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
  if (!usertree) {
    setuptree2(bestree);
    if (njumble > 1)
      setuptree2(bestree2);
  }
  grcategs = (categs > rcategs) ? categs : rcategs;
  allocx(nonodes, grcategs, curtree.nodep, usertree);
  if (!usertree) {
    allocx(nonodes, grcategs, bestree.nodep, 0);
    if (njumble > 1)
      allocx(nonodes, grcategs, bestree2.nodep, 0);
  }
  makevalues2(rcategs, curtree.nodep, endsite, spp, y, alias);
  if (freqsfrom) {
    empiricalfreqs(&freqa, &freqc, &freqg, &freqt, aliasweight, curtree.nodep);
    getbasefreqs(freqa, freqc, freqg, freqt, &freqr, &freqy, &freqar, &freqcy,
                   &freqgr, &freqty, &ttratio, &xi, &xv, &fracchange,
                   freqsfrom, true);
  }
  if (!justwts || firstset)
    fprintf(outfile, "\nTransition/transversion ratio = %10.6f\n\n", ttratio);
}  /* getinput */


void inittable_for_usertree (FILE *intree,char *progname,char *intreename)
{
  /* If there's a user tree, then the ww/zz/wwzz/vvzz elements need
     to be allocated appropriately. */
  long num_comma;
  long i, j;

  /* First, figure out the largest possible furcation, i.e. the number
     of commas plus one */
  countcomma (&intree, &num_comma);
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
  
  tbl = (valrec ***) calloc (1, rcategs * sizeof(valrec **));
  for (i = 0; i < rcategs; i++) {
    tbl[i] = (valrec **) calloc (1, categs*sizeof(valrec *));
    for (j = 0; j < categs; j++)
      tbl[i][j] = (valrec *) calloc (1, sizeof(valrec));
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


void exmake(double lz, long n)
{
  /* pretabulate tables of exponentials so need not do for each site */
  long i;
  double rat;

  for (i = 0; i < categs; i++) {
    rat = rate[i];
    switch (n) {

    case 1:
      expon1i[i] = exp(rat * xi * lz);
      expon1v[i] = exp(rat * xv * lz);
      break;

    case 2:
      expon2i[i] = exp(rat * xi * lz);
      expon2v[i] = exp(rat * xv * lz);
      break;
    }
  }
}  /* exmake */


void alloc_nvd(long num_sibs, nuview_data *local_nvd)
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


void free_nvd(long num_sibs, nuview_data *local_nvd)
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
/* current (modified dnaml) nuview */
{
  long i, j, k, num_sibs, sib_index;
  nuview_data *local_nvd;
  node *sib_ptr, *sib_back_ptr;
  sitelike p_xx;
  double lw;

  /* Figure out how many siblings the current node has */
  num_sibs    = count_sibs (p);
/*xx should this be enabled? */
  /* Recursive calls, should be called for all children */
  sib_ptr = p;
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    if (!(sib_back_ptr == NULL))
      if (!sib_back_ptr->tip && !sib_back_ptr->initialized)
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
    
    if (sib_back_ptr != NULL)
      lw =  -fabs(p->tyme - sib_back_ptr->tyme);
    else
      lw = 0.0;
    

    for (i = 0; i < rcategs; i++)
      for (j = 0; j < categs; j++) {
        tbl[i][j]->ww[sib_index]   = exp(tbl[i][j]->ratxi * lw);
        tbl[i][j]->zz[sib_index]   = exp(tbl[i][j]->ratxv * lw);
        tbl[i][j]->wwzz[sib_index] = tbl[i][j]->ww[sib_index] * tbl[i][j]->zz[sib_index];
        tbl[i][j]->vvzz[sib_index] = (1.0 - tbl[i][j]->ww[sib_index]) *
          tbl[i][j]->zz[sib_index];
      }
  }
  /* edit point: null back pointer protection code */
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
        if (sib_back_ptr != NULL)
          memcpy(local_nvd->xx[sib_index],
               sib_back_ptr->x[i][j],
               sizeof(sitelike));
        else {
          local_nvd->xx[sib_index][0] = 1.0;
          local_nvd->xx[sib_index][(long)C - (long)A] = 1.0;
          local_nvd->xx[sib_index][(long)G - (long)A] = 1.0;
          local_nvd->xx[sib_index][(long)T - (long)A] = 1.0;
        }
      }

      /* Loop 2.2 */
      sib_ptr = p;
      for (sib_index=0; sib_index < num_sibs; sib_index++) {
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


double evaluate(node *p)
{
  contribarr tterm;
  static contribarr like, nulike, clai;
  double sum, sum2, sumc=0, y, lz, y1, z1zz, z1yy, 
		prod12, prod1, prod2, prod3, sumterm, lterm;
  long i, j, k, lai;
  node *q, *r;
  sitelike x1, x2;
  sum = 0.0;

  if (p == curtree.root && (count_sibs(p) == 2)) {
    r = p->next->back;
    q = p->next->next->back;
    y = r->tyme + q->tyme - 2 * p->tyme;
    if (!r->tip && !r->initialized) nuview (r);
    if (!q->tip && !q->initialized) nuview (q);
  } else if (p == curtree.root) {
    /* the next two lines copy tyme and x to p->next.  Normally they are 
       not initialized for an internal node. */
/* assumes bifurcation */
    p->next->tyme = p->tyme;
    nuview(p->next);
    r = p->next;
    q = p->next->back;
    y = fabs(p->next->tyme - q->tyme);
  } else {
    r = p;
    q = p->back;
    if (!r->tip && !r->initialized) nuview (r);
    if (!q->tip && !q->initialized) nuview (q);
    y = fabs(r->tyme - q->tyme);
  }

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
/* orig  p[i][j]    memcpy(x1, p->x[i][j], sizeof(sitelike));*/
/* xx r[i][j]*/      memcpy(x1, r->x[i][j], sizeof(sitelike));
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
    memcpy(contribution[i], clai, sizeof(contribarr));
    if (!auto_ && usertree)
      l0gf[which - 1][i] = lterm;
    sum += aliasweight[i] * lterm;
  }
  if (auto_) {
    for (j = 0; j < rcategs; j++)
      like[j] = 1.0;
    for (i = 0; i < sites; i++) {
      if ((ally[i] > 0) && (location[ally[i]-1] > 0)) {
        sumc = 0.0;
        for (k = 0; k < rcategs; k++)
          sumc += probcat[k] * like[k];
        sumc *= lambda;
        lai = location[ally[i] - 1];
        memcpy(clai, contribution[lai - 1], sizeof(contribarr));
        for (j = 0; j < rcategs; j++)
          nulike[j] = ((1.0 - lambda) * like[j] + sumc) * clai[j];
      } else {
        for (j = 0; j < rcategs; j++)
          nulike[j] = ((1.0 - lambda) * like[j] + sumc);
      }
      memcpy(like, nulike, sizeof(contribarr));
    }
    sum2 = 0.0;
    for (i = 0; i < rcategs; i++)
      sum2 += probcat[i] * like[i];
    sum += log(sum2);
  }
  curtree.likelihood = sum;
  if (auto_ || !usertree)
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


/* original dnamlk nuview */
void oldnuview(node *p)
{
  long i, j;
  double lw, ww1, ww2, zz1, zz2, vv1, vv2, yy1, yy2, sum1, sum2, sumr1, sumr2,
         sumy1, sumy2;
  node *q, *r;
  sitelike xx1, xx2, xx3;

  q = p->next->back;
  r = p->next->next->back;
  if (q != NULL)
    lw = -fabs(p->tyme - q->tyme);
  else
    lw = 0.0;
  exmake(lw, 1L);
  if (r != NULL)
    lw = -fabs(p->tyme - r->tyme);
  else
    lw = 0.0;
  exmake(lw, 2L);
  for (i = 0; i < endsite; i++) {
    for (j = 0; j < categs; j++) {
      ww1 = expon1i[j];
      zz1 = expon1v[j];
      vv1 = 1.0 - ww1;
      yy1 = 1.0 - zz1;
      ww2 = expon2i[j];
      zz2 = expon2v[j];
      vv2 = 1.0 - ww2;
      yy2 = 1.0 - zz2;
      if (q != NULL)
        memcpy(xx1, q->x[i][j], sizeof(sitelike));
      if (r != NULL)
        memcpy(xx2, r->x[i][j], sizeof(sitelike));
      if (q == NULL) {
        xx1[0] = 1.0;
        xx1[(long)C - (long)A] = 1.0;
        xx1[(long)G - (long)A] = 1.0;
        xx1[(long)T - (long)A] = 1.0;
      }
      if (r == NULL) {
        xx2[0] = 1.0;
        xx2[(long)C - (long)A] = 1.0;
        xx2[(long)G - (long)A] = 1.0;
        xx2[(long)T - (long)A] = 1.0;
      }
      sum1 = yy1 * (freqa * xx1[0] + freqc * xx1[(long)C - (long)A] +
            freqg * xx1[(long)G - (long)A] + freqt * xx1[(long)T - (long)A]);
      sum2 = yy2 * (freqa * xx2[0] + freqc * xx2[(long)C - (long)A] +
            freqg * xx2[(long)G - (long)A] + freqt * xx2[(long)T - (long)A]);
      sumr1 = vv1 * (freqar * xx1[0] + freqgr * xx1[(long)G - (long)A]);
      sumr2 = vv2 * (freqar * xx2[0] + freqgr * xx2[(long)G - (long)A]);
      sumy1 = vv1 * (freqcy * xx1[(long)C - (long)A] +
                     freqty * xx1[(long)T - (long)A]);
      sumy2 = vv2 * (freqcy * xx2[(long)C - (long)A] +
                     freqty * xx2[(long)T - (long)A]);
      xx3[0] = (sum1 + zz1 * (ww1 * xx1[0] + sumr1)) *
               (sum2 + zz2 * (ww2 * xx2[0] + sumr2));
      xx3[(long)C - (long)A] =
        (sum1 + zz1 * (ww1 * xx1[(long)C - (long)A] + sumy1)) *
        (sum2 + zz2 * (ww2 * xx2[(long)C - (long)A] + sumy2));
      xx3[(long)G - (long)A] =
        (sum1 + zz1 * (ww1 * xx1[(long)G - (long)A] + sumr1)) *
        (sum2 + zz2 * (ww2 * xx2[(long)G - (long)A] + sumr2));
      xx3[(long)T - (long)A] =
        (sum1 + zz1 * (ww1 * xx1[(long)T - (long)A] + sumy1)) *
        (sum2 + zz2 * (ww2 * xx2[(long)T - (long)A] + sumy2));
      memcpy(p->x[i][j], xx3, sizeof(sitelike));
    }
  }
}  /* oldnuview */


void getthree(node *p)
{
  /* compute likelihood, slope, curvature at a new triple of points */
  double tt, td;
  node *sdown, *sib_ptr, *sib_back_ptr;
  long num_sibs, i;

  sdown = curtree.nodep[p->index - 1];
  sib_ptr = sdown;
  sdown = sdown->back;
  td = fabs(tdelta);
  tt = p->tyme;

  num_sibs = count_sibs(p);
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    if (sib_back_ptr->tyme - tt < td + epsilon)
      td = (sib_back_ptr->tyme - tt) / 2.0;
  }

  if (sdown != NULL) {
    if (tt - sdown->tyme < td + epsilon)
      td = (tt - sdown->tyme) / 2.0;
  }
  if (-tt > epsilon && td < epsilon)
    td = epsilon;
  if (tt == -td)
    td += epsilon;
  p->tyme = tt + td;
  x[0] = tt + td;
  nuview(p);
  lnl[0] = evaluate(p);
  p->tyme = tt - td;
  x[2] = tt - td;
  nuview(p);
  lnl[2] = evaluate(p);
  p->tyme = tt;
  x[1] = tt;
  nuview(p);
  lnl[1] = evaluate(p);
}  /* getthree */


void makenewv(node *p)
{
  /* improve a node time */
  long it, imin, imax, i, num_sibs;
  double tt, tfactor, tlow, thigh, oldlike, ymin, ymax, s32, s21, yold;
  boolean done, already;
  node *s, *sdown, *sib_ptr, *sib_back_ptr;

  s = curtree.nodep[p->index - 1];
  sdown = s->back;
  if (s == curtree.root)
    tlow = -10.0;
  else
    tlow = sdown->tyme;

  sib_ptr = s;
  num_sibs = count_sibs(p);

  thigh = s->next->back->tyme;
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
      if (sib_back_ptr->tyme < thigh)
        thigh = sib_back_ptr->tyme;
  }
  done = (thigh - tlow < 2.0*epsilon);
  it = 1;
  if (s != curtree.root)
    tdelta = (thigh - tlow) / 10.0;
  else
    tdelta = (thigh - s->tyme) / 5.0;
  tfactor = 1.0;
  if (!done)
    getthree(s);
  while (it < iterations && !done) {
    tt = s->tyme;
    ymax = lnl[0];
    imax = 1;
    for (i = 2; i <= 3; i++) {
      if (lnl[i - 1] > ymax) {
        ymax = lnl[i - 1];
        imax = i;
      }
    }
    if (imax != 2) {
      ymax = x[1];
      x[1] = x[imax - 1];
      x[imax - 1] = ymax;
      ymax = lnl[1];
      lnl[1] = lnl[imax - 1];
      lnl[imax - 1] = ymax;
    }
    tt = x[1];
    oldlike = lnl[1];
    yold = tt;
    s32 = (lnl[2] - lnl[1]) / (x[2] - x[1]);
    s21 = (lnl[1] - lnl[0]) / (x[1] - x[0]);
    if (fabs(x[2] - x[0]) > epsilon)
      curv = (s32 - s21) / ((x[2] - x[0]) / 2);
    else
      curv = 0.0;
    slope = (s32 + s21) / 2 - curv *          (x[2] - 2 * x[1] + x[0]) / 4;
    if (curv >= 0.0) {
      if (slope < 0)
        tdelta = -fabs(tdelta);
      else
        tdelta = fabs(tdelta);
    } else
      tdelta = -(tfactor * slope / curv);
    if (tt + tdelta <= tlow + epsilon)
      tdelta = (tlow - tt) / 2;
    if (tt + tdelta >= thigh - epsilon)
      tdelta = (thigh - tt) / 2;
    tt += tdelta;
    done = (fabs(yold - tt) < epsilon || fabs(tdelta) < epsilon);
    s->tyme = tt;
    nuview(s);
    lnlike = evaluate(s);
    ymin = lnl[0];
    imin = 1;
    for (i = 2; i <= 3; i++) {
      if (lnl[i - 1] < ymin) {
        ymin = lnl[i - 1];
        imin = i;
      }
    }
    already = (tt == x[0]) || (tt == x[1]) || (tt == x[2]);
    if (!already && ymin < lnlike) {
      x[imin - 1] = tt;
      lnl[imin - 1] = lnlike;
    }
    if (already || lnlike < oldlike) {
      tt = x[1];
      tfactor /= 2;
      tdelta /= 2;
      curtree.likelihood = oldlike;
      lnlike = oldlike;
    } else
      tfactor = 1.0;


    sib_ptr = p;
    num_sibs = count_sibs(p);
    p->tyme = tt;
    for (i=0 ; i < num_sibs; i++) {
      sib_ptr      = sib_ptr->next;
      sib_ptr->tyme = tt;
    }

    sib_ptr = p;
    nuview(p);
    for (i=0 ; i < num_sibs; i++) {
      sib_ptr      = sib_ptr->next;
      nuview(sib_ptr);
    }
    
    it++;
  }
  sib_ptr = p;
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    inittrav (sib_ptr);
  }
  smoothed = smoothed && done;
}  /* makenewv */


void update(node *p)
{
  node *sib_ptr, *sib_back_ptr;
  long i, num_sibs;
  
  /* improve time and recompute views at a node */
  if (p == NULL)
    return;
  if (p->back != NULL) {
    if (!p->back->tip && !p->back->initialized)
      nuview(p->back);
  }

  sib_ptr = p;
  num_sibs = count_sibs(p);
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    if (sib_back_ptr != NULL) {
      if (!sib_back_ptr->tip && !sib_back_ptr->initialized)
        nuview(sib_back_ptr);
    }
  }

  if ((!usertree) || (usertree && !lengths) || p->iter) {
    makenewv(p);
    return;
  }
  nuview(p);

  sib_ptr = p;
  num_sibs = count_sibs(p);
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    nuview(sib_ptr);  
  }
}  /* update */


void smooth(node *p)
{
  node *sib_ptr;
  long i, num_sibs;

  if (p == NULL)
    return;
  if (p->tip)
    return;

  update(p);

  smoothed = false;
  sib_ptr = p;
  num_sibs = count_sibs(p);
  for (i=0; i < num_sibs; i++) {
    sib_ptr = sib_ptr->next;
    if (polishing || (smoothit && !smoothed)) {
      smooth(sib_ptr->back);
      p->initialized = false;
      sib_ptr->initialized = false;
    }
    update(p);  
  }
}  /* smooth */


void restoradd(node *below, node *newtip, node *newfork, double prevtyme)
{
/* restore "new" tip and fork to place "below".  restore tymes */
/* assumes bifurcation */
  hookup(newfork, below->back);
  hookup(newfork->next, below);
  hookup(newtip, newfork->next->next);
  curtree.nodep[newfork->index-1] = newfork;
  newfork->tyme = prevtyme;
/* assumes bifurcations */
  newfork->next->tyme = prevtyme;
  newfork->next->next->tyme = prevtyme;
} /* restoradd */


void dnamlk_add(node *below, node *newtip, node *newfork, boolean tempadd)
{
  /* inserts the nodes newfork and its descendant, newtip, into the tree. */
  long i;
  boolean done;
  node *p;

  below = curtree.nodep[below->index - 1];
  newfork = curtree.nodep[newfork->index-1];
  newtip = curtree.nodep[newtip->index-1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  below->back = newfork->next->next;
  newfork->next->next->back = below;
  newfork->next->back = newtip;
  newtip->back = newfork->next;
  if (newtip->tyme < below->tyme)
    p = newtip;
  else p = below;
  newfork->tyme = p->tyme;
  if (curtree.root == below)
    curtree.root = newfork;
  if (newfork->back != NULL) {
    if (p->tyme > newfork->back->tyme)
      newfork->tyme = (p->tyme + newfork->back->tyme) / 2.0;
    else newfork->tyme = p->tyme - epsilon;
    newfork->next->tyme = newfork->tyme;
    newfork->next->next->tyme = newfork->tyme;
    do {
      p = curtree.nodep[p->back->index - 1];
      done = (p == curtree.root);
      if (!done)
        done = (curtree.nodep[p->back->index - 1]->tyme < p->tyme - epsilon);
      if (!done) {
        curtree.nodep[p->back->index - 1]->tyme = p->tyme - epsilon;
        curtree.nodep[p->back->index - 1]->next->tyme = p->tyme - epsilon;
        curtree.nodep[p->back->index - 1]->next->next->tyme = p->tyme - epsilon;
      }
    } while (!done);
  } else {
      newfork->tyme = newfork->tyme - 2*epsilon;
      newfork->next->tyme = newfork->tyme;
      newfork->next->next->tyme = newfork->tyme;
    }
  inittrav(newtip);
  inittrav(newtip->back);
  smoothed = false;
  i = 1;
  while (i < smoothings && !smoothed) {
    smoothed = true;
    smooth(newfork);
    smooth(newfork->back);
    i++;
  }
}  /* dnamlk_add */


void dnamlk_re_move(node **item, node **fork, boolean tempadd)
{
  /* removes nodes item and its ancestor, fork, from the tree.
    the new descendant of fork's ancestor is made to be
    fork's second descendant (other than item).  Also
    returns pointers to the deleted nodes, item and fork */
  node *p, *q;
  long i;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *item = curtree.nodep[(*item)->index-1];
  *fork = curtree.nodep[(*item)->back->index - 1];
  if (curtree.root == *fork) {
    if (*item == (*fork)->next->back)
      curtree.root = (*fork)->next->next->back;
    else
      curtree.root = (*fork)->next->back;
  }
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
/* debug replace by hookup calls?  Does that have NULL protection? */
  if (p != NULL)
    p->back = q;
  if (q != NULL)
    q->back = p;
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
  inittrav(p);
  inittrav(q);
  if (tempadd)
    return;
  i = 1;
  while (i <= smoothings) {
    smooth(q);
    if (smoothit)
      smooth(q->back);
    i++;
  }
}  /* dnamlk_re_move */


void tryadd(node *p, node **item, node **nufork)
{  /* temporarily adds one fork and one tip to the tree.
    if the location where they are added yields greater
    likelihood than other locations tested up to that
    time, then keeps that location as there */

  long grcategs;
  grcategs = (categs > rcategs) ? categs : rcategs;
  
  dnamlk_add(p, *item, *nufork, true);
  like = evaluate(p);
  if (lastsp) {
      if (like >= bestyet) 
            copy_(&curtree, &bestree, nonodes, grcategs);
  }
  if (like > bestyet) {
    bestyet = like;
    there = p;
  }
  dnamlk_re_move(item, nufork, true);
}  /* tryadd */


void addpreorder(node *p, node *item_, node *nufork_, boolean contin,
			boolean continagain)
{
  /* traverses a binary tree, calling function tryadd
    at a node before calling tryadd at its descendants */
  node *item, *nufork;

  item = item_;
  nufork = nufork_;
  if (p == NULL)
    return;
  tryadd(p, &item, &nufork);
  contin = continagain;
  if ((!p->tip) && contin) {
    addpreorder(p->next->back, item, nufork, contin, continagain);
    addpreorder(p->next->next->back, item, nufork, contin, continagain);
  }
}  /* addpreorder */


void tryrearr(node *p, boolean *success)
{
  /* evaluates one rearrangement of the tree.
    if the new tree has greater likelihood than the old
    one sets success = TRUE and keeps the new tree.
    otherwise, restores the old tree */
  node *frombelow, *whereto, *forknode;
  double oldlike, prevtyme;
  boolean wasonleft;

  if (p == curtree.root)
    return;
  forknode = curtree.nodep[p->back->index - 1];
  if (forknode == curtree.root)
    return;
  oldlike = bestyet;
  prevtyme = forknode->tyme;
/* the following statement presumes bifurcating tree */
  if (forknode->next->back == p) {
    frombelow = forknode->next->next->back;
    wasonleft = true;
  }
  else {
    frombelow = forknode->next->back;
    wasonleft = false;
  }
  whereto = curtree.nodep[forknode->back->index - 1];
  dnamlk_re_move(&p, &forknode, true);
  dnamlk_add(whereto, p, forknode, true);
  like = evaluate(p);
  if (like <= oldlike) {
    dnamlk_re_move(&p, &forknode, true);
    restoradd(frombelow, p, forknode, prevtyme);
    if (wasonleft && (forknode->next->next->back == p)) {
       hookup (forknode->next->back, forknode->next->next);
       hookup (forknode->next, p);
    }
    curtree.likelihood = oldlike;
    inittrav(forknode);
    inittrav(forknode->next);
    inittrav(forknode->next->next);
  } else {
    (*success) = true;
    bestyet = like;
  }
}  /* tryrearr */


void repreorder(node *p, boolean *success)
{
  /* traverses a binary tree, calling function tryrearr
    at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p, success);
  if (p->tip)
    return;
  if (!(*success))
    repreorder(p->next->back, success);
  if (!(*success))
    repreorder(p->next->next->back, success);
}  /* repreorder */


void rearrange(node **r)
{
  /* traverses the tree (preorder), finding any local
    rearrangement which increases the likelihood.
    if traversal succeeds in increasing the tree's
    likelihood, function rearrange runs traversal again */
  boolean success;
  success = true;
  while (success) {
    success = false;
    repreorder(*r, &success);
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
  case hslength:
    break;
  case hsnolength:
    break;
  case treewt:
    break;
  case unittrwt:
    break;
  }
} /* initdnamlnode */


void tymetrav(node *p, double *x)
{
  /* set up times of nodes */
  node *sib_ptr;
  long i, num_sibs;

  if (!p->tip) {
    sib_ptr  = p;
    num_sibs = count_sibs(p);
    for (i=0; i < num_sibs; i++) {
      sib_ptr = sib_ptr->next;
      tymetrav(sib_ptr->back, x);
    }
  } else
    (*x)     = 0.0;
    p->tyme  = (*x);
    (*x)    -= p->v;
}  /* tymetrav */


void dnamlk_coordinates(node *p, long *tipy)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last, *pp1 =NULL, *pp2 =NULL;
  long num_sibs, p1, p2, i;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = (*tipy);
    p->ymin   = (*tipy);
    p->ymax   = (*tipy);
    (*tipy)  += down;
    return;
  }
  q = p->next;
  do {
    dnamlk_coordinates(q->back, tipy);
    q = q->next;
  } while (p != q);
  num_sibs = count_sibs(p);
  p1 = (long)((num_sibs+1)/2.0);
  p2 = (long)((num_sibs+2)/2.0);
  i = 1;
  q = p->next;
  first  = q->back;
  do {
    if (i == p1) pp1 = q->back;
    if (i == p2) pp2 = q->back;
    last = q->back;
    q = q->next;
    i++;
  } while (q != p);
  p->xcoord = (long)(0.5 - over * p->tyme);
  p->ycoord = (pp1->ycoord + pp2->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* dnamlk_coordinates */


void dnamlk_drawline(long i, double scale)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first =NULL, *last =NULL;
  long n, j;
  boolean extra, done;

  p = curtree.root;
  q = curtree.root;
  extra = false;
  if ((long)(p->ycoord) == i) {
    if (p->index - spp >= 10)
      fprintf(outfile, "-%2ld", p->index - spp);
    else
      fprintf(outfile, "--%ld", p->index - spp);
    extra = true;
  } else
    fprintf(outfile, "  ");
  do {
    if (!p->tip) {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || r == p));
      first = p->next->back;
      r = p->next;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p == q);
    n = (long)(scale * ((long)(p->xcoord) - (long)(q->xcoord)) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)(q->ycoord) == i && !done) {
      if (p->ycoord != q->ycoord)
        putc('+', outfile);
      else
        putc('-', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip) {
      if ((long)(last->ycoord) > i && (long)(first->ycoord) < i && 
           i != (long)(p->ycoord)) {
        putc('!', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
    }
    if (p != q)
      p = q;
  } while (!done);
  if ((long)(p->ycoord) == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* dnamlk_drawline */


void dnamlk_printree()
{
 /* prints out diagram of the tree */
  long tipy;
  double scale;
  long i;
  node *p;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  dnamlk_coordinates(curtree.root, &tipy);
  p = curtree.root;
  while (!p->tip)
    p = p->next->back;
  scale = 1.0 / (long)(p->tyme - curtree.root->tyme + 1.000);
  putc('\n', outfile);
  for (i = 1; i <= tipy - down; i++)
    dnamlk_drawline(i, scale);
  putc('\n', outfile);
}  /* dnamlk_printree */


void describe(node *p)
{
  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;
  double v;

  if (p == curtree.root)
    fprintf(outfile, " root         ");
  else
    fprintf(outfile, "%4ld          ", p->back->index - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", p->index - spp);
  if (p != curtree.root) {
    fprintf(outfile, "%11.5f", fracchange * (p->tyme - curtree.root->tyme));
    v = fracchange * (p->tyme - curtree.nodep[p->back->index - 1]->tyme);
    fprintf(outfile, "%13.5f", v);
  }
  putc('\n', outfile);
  if (!p->tip) {

    sib_ptr = p;
    num_sibs = count_sibs(p);
    for (i=0 ; i < num_sibs; i++) {
      sib_ptr      = sib_ptr->next;
      sib_back_ptr = sib_ptr->back;
      describe(sib_back_ptr);
    }
  }
}  /* describe */


void oldsummarize()
{
  long i, j, mx;
  double mode, sum;
  double *like, *nulike;
  long **mp;

  like    = (double *)Malloc(categs * sizeof(double));
  nulike  = (double *)Malloc(categs * sizeof(double));
  mp      = (long **)Malloc(sites * sizeof(long *));
  for (i=0;i<=sites-1;++i)
     mp[i] = (long *)Malloc(sizeof(long)*categs);
  fprintf(outfile, "\nLn Likelihood = %11.5f\n\n", curtree.likelihood);
  fprintf(outfile, " Ancestor      Node      Node Height     Length\n");
  fprintf(outfile, " --------      ----      ---- ------     ------\n");
  describe(curtree.root);
  putc('\n', outfile);
  if (ctgry && categs > 1) {
    for (i = 0; i < categs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--) {
      sum = 0.0;
      for (j = 0; j < categs; j++) {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
        mp[i][j] = j + 1;
        for (k = 1; k <= categs; k++) {
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
      for (j = 0; j < categs; j++)
        nulike[j] /= sum;
      memcpy(like, nulike, categs * sizeof(double));
    }
    mode = 0.0;
    mx = 1;
    for (i = 1; i <= categs; i++) {
      if (probcat[i - 1] * like[i - 1] > mode) {
        mx = i;
        mode = probcat[i - 1] * like[i - 1];
      }
    }
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
  }
  putc('\n', outfile);
  putc('\n', outfile);
  free(like);
  free(nulike);
  for (i=0;i<sites;++i)
    free(mp[i]);
  free(mp);
}  /* oldsummarize */


void summarize()
{
  long i, j, mx, mx0;
  double mode, sum;
  double *like, *nulike;
  double **marginal;
  long **mp;

  like    = (double *)Malloc(categs * sizeof(double));
  nulike  = (double *)Malloc(categs * sizeof(double));
  mp      = (long **)Malloc(sites * sizeof(long *));
  for (i=0;i<=sites-1;++i)
     mp[i] = (long *)Malloc(sizeof(long)*categs);
  fprintf(outfile, "\nLn Likelihood = %11.5f\n\n", curtree.likelihood);
  fprintf(outfile, " Ancestor      Node      Node Height     Length\n");
  fprintf(outfile, " --------      ----      ---- ------     ------\n");
  describe(curtree.root);
  putc('\n', outfile);
/**/
  if (rctgry && rcategs > 1) {
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++) {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
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
    marginal = (double **) calloc (1, sites*sizeof(double *));
    for (i = 0; i < sites; i++)
      marginal[i] = (double *) calloc (1, rcategs*sizeof(double));
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++) {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
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
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
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
/**/
  putc('\n', outfile);
  putc('\n', outfile);
  free(like);
  free(nulike);
  for (i=0;i<sites;++i)
    free(mp[i]);
  free(mp);
}  /* summarize */


void dnamlk_treeout(node *p)
{
  /* write out file with representation of final tree */
  node *sib_ptr;
  long i, n, w, num_sibs;
  Char c;
  double x;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    col += n;
  } else {
    sib_ptr = p;
    num_sibs = count_sibs(p);
    putc('(', outtree);
    col++;

    for (i=0; i < (num_sibs - 1); i++) {
      sib_ptr = sib_ptr->next;
      dnamlk_treeout(sib_ptr->back);
      putc(',', outtree);
      col++;
      if (col > 55) {
        putc('\n', outtree);
        col = 0;
      }
    }
    sib_ptr = sib_ptr->next;
    dnamlk_treeout(sib_ptr->back);
    putc(')', outtree);
    col++;
  }
  if (p == curtree.root) {
    fprintf(outtree, ";\n");
    return;
  }
  x = fracchange * (p->tyme - curtree.nodep[p->back->index - 1]->tyme);
  if (x > 0.0)
    w = (long)(0.4342944822 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.4342944822 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  fprintf(outtree, ":%*.5f", (int)(w + 7), x);
  col += w + 8;
}  /* dnamlk_treeout */


void nodeinit(node *p)
{
  /* set up times at one node */
  node *sib_ptr, *sib_back_ptr;
  long i, num_sibs;
  double lowertyme;

  sib_ptr = p;
  num_sibs = count_sibs(p);

  /* lowertyme = lowest of children's times */
  lowertyme = p->next->back->tyme;
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    if (sib_back_ptr->tyme < lowertyme)
      lowertyme = sib_back_ptr->tyme;
  }

  p->tyme = lowertyme - 0.1;

  sib_ptr = p;
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    sib_ptr->tyme = p->tyme;
    sib_back_ptr->v = sib_back_ptr->tyme - p->tyme;
    sib_ptr->v = sib_back_ptr->v;
  }
}  /* nodeinit */


void oldnodeinit(node *p)
{
  /* set up times at one node */
  node *q, *r;
  double lowertyme;

  q = p->next->back;
  r = p->next->next->back;

/* lowertyme = lowest of children's times */
  lowertyme = q->tyme;
  if (r->tyme < q->tyme)
    lowertyme = r->tyme;

/* keep */
  p->tyme = lowertyme - 0.1;

/* Set tyme of node ring members to that of p */
  p->next->tyme = p->tyme;
  p->next->next->tyme = p->tyme;

/* for all (in this case, both) children:
    child->v = child->tyme - p->tyme
    (ring member whose back pointer points at child)->v = child->v
*/
  q->v = q->tyme - p->tyme;
  p->next->v = q->v;
  r->v = r->tyme - p->tyme;
  p->next->next->v = r->v;
}  /* oldnodeinit */


void initrav(node *p)
{

  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;

  /* traverse to set up times throughout tree */
  if (p->tip)
    return;

  sib_ptr = p;
  num_sibs = count_sibs(p);
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    initrav(sib_back_ptr);
  }

  nodeinit(p);
}  /* initrav */


void travinit(node *p)
{
  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;

  /* traverse to set up initial values */
  if (p == NULL)
    return;
  if (p->tip)
    return;
  if (p->initialized)
    return;


  sib_ptr = p;
  num_sibs = count_sibs(p);
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    travinit(sib_back_ptr);
  }

  nuview(p);
  p->initialized = true;
}  /* travinit */


void travsp(node *p)
{
  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;

  /* traverse to find tips */
  if (p == curtree.root)
    travinit(p);
  if (p->tip)
    travinit(p->back);
  else {
    sib_ptr = p;
    num_sibs = count_sibs(p);
    for (i=0 ; i < num_sibs; i++) {
      sib_ptr      = sib_ptr->next;
      sib_back_ptr = sib_ptr->back;
      travsp(sib_back_ptr);
    }
  }
}  /* travsp */


void treevaluate()
{
  /* evaluate likelihood of tree, after iterating branch lengths */
  long i, j,  num_sibs;
  node *sib_ptr, *sib_back_ptr;
  double dummy;

  polishing = true;
  smoothit = true;
  for (i = 0; i < spp; i++)
    curtree.nodep[i]->initialized = false;
  for (i = spp; i < nonodes; i++) {
    sib_ptr = curtree.nodep[i];
    sib_ptr->initialized = false;
    num_sibs = count_sibs(sib_ptr);
    for (j=0 ; j < num_sibs; j++) {
      sib_ptr      = sib_ptr->next;
      sib_back_ptr = sib_ptr->back;
      sib_ptr->initialized = false;
    }

  }
  if (!lengths)
    initrav(curtree.root);
  travsp(curtree.root);
  for (i = 1; i <= smoothings * 4; i++)
    smooth(curtree.root);
  dummy = evaluate(curtree.root);
}  /* treevaluate */


void maketree()
{
  /* constructs a binary tree from the pointers in curtree.nodep,
     adds each node at location which yields highest likelihood
     then rearranges the tree for greatest likelihood */

  long i, j, numtrees;
  double bestlike, gotlike, x;
  node *item, *nufork, *dummy, *q, *root=NULL;
  boolean dummy_haslengths, dummy_first, goteof;
  long nextnode;
  long grcategs;
  pointarray dummy_treenode;

  grcategs = (categs > rcategs) ? categs : rcategs;

  inittable();  

  if (!usertree) {
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    curtree.root = curtree.nodep[spp];
    curtree.root->back = NULL;
    for (i = 0; i < spp; i++)
       curtree.nodep[i]->back = NULL;
    for (i = spp; i < nonodes; i++) {
       q = curtree.nodep[i];
       q->back = NULL;
       while ((q = q->next) != curtree.nodep[i])
         q->back = NULL;
    }
    polishing = false;
    dnamlk_add(curtree.nodep[enterorder[0] - 1], curtree.nodep[enterorder[1] - 1],
        curtree.nodep[spp], false);
    if (progress) {
      printf("\nAdding species:\n");
      writename(0, 2, enterorder);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    lastsp = false;
    smoothit = false;
    for (i = 3; i <= spp; i++) {
      bestree.likelihood = spp*sites*(freqa*log(freqa)+freqg*log(freqg) +
                                freqc*log(freqc)+freqt*log(freqt));
      bestyet = spp*sites*(freqa*log(freqa)+freqg*log(freqg) +
                                freqc*log(freqc)+freqt*log(freqt));
      there = curtree.root;
      item = curtree.nodep[enterorder[i - 1] - 1];
      nufork = curtree.nodep[spp + i - 2];
      lastsp = (i == spp);
      addpreorder(curtree.root, item, nufork, true, true);
      dnamlk_add(there, item, nufork, false);
      like = evaluate(curtree.root);
      copy_(&curtree, &bestree, nonodes, grcategs);
      rearrange(&curtree.root);
      if (curtree.likelihood > bestree.likelihood) {
        copy_(&curtree, &bestree, nonodes, grcategs);
      }
      if (progress) {
        writename(i - 1, 1, enterorder);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
      if (lastsp && global) {
        if (progress) {
          printf("Doing global rearrangements\n");
          printf("  !");
          for (j = 1; j <= nonodes; j++)
            putchar('-');
          printf("!\n");
        }
        bestlike = bestyet;
        do {
          if (progress)
            printf("   ");
          gotlike = bestlike;
          for (j = 0; j < nonodes; j++) {
            bestyet = spp*sites*(freqa*log(freqa)+freqg*log(freqg) +
                                freqc*log(freqc)+freqt*log(freqt));
            item = curtree.nodep[j];
            if (item != curtree.root) {
              nufork = curtree.nodep[curtree.nodep[j]->back->index - 1];
              dnamlk_re_move(&item, &nufork, false);
              there = curtree.root;
              addpreorder(curtree.root, item, nufork, true, true);
              dnamlk_add(there, item, nufork, false);
            }
            if (progress) {
              putchar('.');
              fflush(stdout);
            }
          }
          if (progress)
            putchar('\n');
        } while (bestlike < gotlike);
      }
    }
    if (njumble > 1 && lastsp) {
      for (i = 0; i < spp; i++ )
        dnamlk_re_move(&curtree.nodep[i], &dummy, false);
      if (jumb == 1 || bestree2.likelihood < bestree.likelihood)
        copy_(&bestree, &bestree2, nonodes, grcategs);
    }
    if (jumb == njumble) {
      if (njumble > 1)
        copy_(&bestree2, &curtree, nonodes, grcategs);
      else copy_(&bestree, &curtree, nonodes, grcategs);
      fprintf(outfile, "\n\n");
      treevaluate();
      curtree.likelihood = evaluate(curtree.root);
      dnamlk_printree();
      summarize();
      if (trout) {
        col = 0;
        dnamlk_treeout(curtree.root);
      }
    } 
  } else {
    openfile(&intree, INTREE, "input tree file", "r", progname, intreename);
    inittable_for_usertree (intree, progname, intreename);
    numtrees = countsemic(&intree);
    l0gl = (double *)Malloc(numtrees * sizeof(double));
    l0gf = (double **)Malloc(numtrees * sizeof(double *));
    for (i=0;i<numtrees;++i)
      l0gf[i] = (double *)Malloc(endsite * sizeof(double));
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    fprintf(outfile, "\n\n");
    which = 1;
    while (which <= numtrees) {
      
      /* These initializations required each time through the loop
         since multiple trees require re-initialization */
      dummy_haslengths = true;
      nextnode         = 0;
      dummy_first      = true;
      goteof           = false;

      treeread(intree, &root, dummy_treenode, &goteof, &dummy_first,
               curtree.nodep, INTREE, "Dnamlk", &nextnode,
               &dummy_haslengths, &grbg, initdnamlnode);

      nonodes = nextnode;
      
      root = curtree.nodep[root->index - 1];
      curtree.root = root;

      if (lengths)
        tymetrav(curtree.root, &x);

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
      dnamlk_printree();
      summarize();
      if (trout) {
        col = 0;
        dnamlk_treeout(curtree.root);
      }
      which++;
    }      

    FClose(intree);
    if (!auto_ && numtrees > 1 && weightsum > 1 )
      standev2(numtrees, weightsum, maxwhich, 0, endsite, maxlogl, l0gl, l0gf,
               aliasweight);
  }
  
  if (jumb == njumble) {
    if (progress) {
      printf("\nOutput written to output file\n\n");
      if (trout)
        printf("Tree also written onto file\n\n");
    }
    free(contribution);
    freex(nonodes, curtree.nodep);
    if (!usertree) {
      freex(nonodes, bestree.nodep);
      if (njumble > 1)
        freex(nonodes, bestree2.nodep);
    }
  }
  free(root);
} /* maketree */

/*?? Dnaml has a clean-up function for freeing memory, closing files, etc.
     Put one here too? */

int main(int argc, Char *argv[])
{  /* DNA Maximum Likelihood with molecular clock */

#ifdef MAC
  argc = 1;		/* macsetup("Dnamlk", "Dnamlk");	*/
  argv[0] = "Dnamlk";
#endif
#ifdef WIN32
  phySetConsoleAttributes();
  phyClearScreen();
#endif
  strcpy(progname, argv[0]);
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  datasets = 1;
  mulsets = false;
  firstset = true;
  doinit();

  ttratio0    = ttratio;
  if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);
  if (ctgry)
    openfile(&catfile, CATFILE, "categories file", "r", argv[0], catfilename);
  if (weights || justwts)
   openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);
  for (ith = 1; ith <= datasets; ith++) {
    ttratio = ttratio0;
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if (progress)
        printf("\nData set # %ld:\n", ith);
    }
    getinput();

    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
  }
  FClose(infile);  
  FClose(outfile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* DNA Maximum Likelihood with molecular clock */

