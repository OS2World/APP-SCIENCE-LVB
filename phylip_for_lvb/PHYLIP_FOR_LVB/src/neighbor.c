
#include "phylip.h"
#include "dist.h"

/* version 3.6. (c) Copyright 1993-2000 by the University of Washington.
   Written by Mary Kuhner, Jon Yamato, Joseph Felsenstein, Akiko Fuseki,
   Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

extern void alloctree() ;
extern void treeoutr() ;

#ifndef OLDC
/* function prototypes */
void getoptions(void);
void allocrest(void);
void doinit(void);
void inputoptions(void);
void getinput(void);
void describe(node *);
void summarize(void);
void nodelabel(boolean);
void jointree(void);
void maketree(void);
/* function prototypes */
#endif


extern FILE *infile, *outfile, *outtree;
extern long spp;
extern naym *nayme;

Char infilename[100], outfilename[100], outtreename[100];
long nonodes2, outgrno, col, datasets, ith;
long inseed;
vector *x;
intvector *reps;
boolean jumble, lower, upper, outgropt, replicates, trout,
               printdata, progress, treeprint, mulsets, njoin;
tree curtree;
longer seed;
long *enterorder;
Char progname[20];

/* variables for maketree, propagated globally for C version: */
node **cluster;


void getoptions()
{
  /* interactively set options */
  long i;
  long inseed0 = 0;
  Char ch;

  fprintf(outfile, "\nNeighbor-Joining/UPGMA method version %s\n\n",VERSION);
  putchar('\n');
  jumble = false;
  lower = false;
  outgrno = 1;
  outgropt = false;
  replicates = false;
  trout = true;
  upper = false;
  printdata = false;
  progress = true;
  treeprint = true;
  njoin = true;
  for(;;) {
    printf(ansi ? "\033[2J\033[H" : "\n");
    printf("\nNeighbor-Joining/UPGMA method version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  N       Neighbor-joining or UPGMA tree?  %s\n",
           (njoin ? "Neighbor-joining" : "UPGMA"));
    if (njoin) {
      printf("  O                        Outgroup root?");
      if (outgropt)
        printf("  Yes, at species number%3ld\n", outgrno);
      else
        printf("  No, use as outgroup species%3ld\n", outgrno);
    }
    printf("  L         Lower-triangular data matrix?  %s\n",
           (lower ? "Yes" : "No"));
    printf("  R         Upper-triangular data matrix?  %s\n",
           (upper ? "Yes" : "No"));
    printf("  S                        Subreplicates?  %s\n",
           (replicates ? "Yes" : "No"));
    printf("  J     Randomize input order of species?");
    if (jumble)
      printf("  Yes (random number seed =%8ld)\n", inseed0);
    else
      printf("  No. Use input order\n");
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           (ibmpc ? "IBM PC" : ansi ? "ANSI" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    printf("\n\n  Y to accept these or type the letter for one to change\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if  (ch == 'Y')
      break;
    if (strchr("NJOULRSM01234",ch) != NULL){
      switch (ch) {
	
      case 'J':
	jumble = !jumble;
 	if (jumble) {
 	  do {
 	    printf("Random number seed (must be odd)?\n");
#ifdef WIN32
            phyFillScreenColor();
#endif
 	    scanf("%ld%*[^\n]", &inseed);
 	    getchar();
 	  } while (!(inseed & 1));
 	  inseed0 = inseed;
 	  for (i = 0; i <= 5; i++)
 	    seed[i] = 0;
 	  i = 0;
 	  do {
 	    seed[i] = inseed & 63;
 	    inseed /= 64;
 	    i++;
 	  } while (inseed != 0);
 	}
	break;
	
      case 'L':
	lower = !lower;
	break;
	
      case 'O':
	outgropt = !outgropt;
	if (outgropt)
	  initoutgroup(&outgrno, spp);
        else
          outgrno = 1;
	break;
	
      case 'R':
	upper = !upper;
	break;
	
      case 'S':
	replicates = !replicates;
	break;
	
      case 'N':
	njoin = !njoin;
	break;
	
      case 'M':
	mulsets = !mulsets;
	if (mulsets)
	  initdatasets(&datasets);
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
}  /* getoptions */


void allocrest()
{
  long i;

  x = (vector *)Malloc(spp*sizeof(vector));
  for (i = 0; i < spp; i++)
    x[i] = (vector)Malloc(spp*sizeof(double));
  reps = (intvector *)Malloc(spp*sizeof(intvector));
  for (i = 0; i < spp; i++)
    reps[i] = (intvector)Malloc(spp*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  enterorder = (long *)Malloc(spp*sizeof(long));
  cluster = (node **)Malloc(spp*sizeof(node *));
}  /* allocrest */


void doinit()
{
  /* initializes variables */
  node *p;

  inputnumbers2(&spp, &nonodes2, 2);
  nonodes2 += (njoin ? 0 : 1);
  getoptions();
  alloctree(&curtree.nodep, nonodes2+1);
  p = curtree.nodep[nonodes2]->next->next;
  curtree.nodep[nonodes2]->next = curtree.nodep[nonodes2];
  free(p);
  allocrest();
}  /* doinit */


void inputoptions()
{
  /* read options information */

  if (ith != 1)
    samenumsp2(ith);
  putc('\n', outfile);
  if (njoin)
    fprintf(outfile, " Neighbor-joining method\n");
  else
    fprintf(outfile, " UPGMA method\n");
  fprintf(outfile, "\n Negative branch lengths allowed\n\n");
}  /* inputoptions */


void getinput()
{
  /* reads the input data */
  inputoptions();
}  /* getinput */


void describe(node *p)
{
  /* print out information for one branch */
  long i;
  node *q;

  q = p->back;
  fprintf(outfile, "%4ld          ", q->index - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "%15.5f\n", q->v);
  if (!p->tip) {
    describe(p->next->back);
    describe(p->next->next->back);
  }
}  /* describe */


void summarize()
{
  /* print out branch lengths etc. */
  putc('\n', outfile);
  if (njoin) {
    fprintf(outfile, "remember:");
    if (outgropt)
      fprintf(outfile, " (although rooted by outgroup)");
    fprintf(outfile, " this is an unrooted tree!\n");
  }
  fprintf(outfile, "\nBetween        And            Length\n");
  fprintf(outfile, "-------        ---            ------\n");
  describe(curtree.start->next->back);
  describe(curtree.start->next->next->back);
  if (njoin)
    describe(curtree.start->back);
  fprintf(outfile, "\n\n");
}  /* summarize */


void nodelabel(boolean isnode)
{
  if (isnode)
    printf("node");
  else
    printf("OTU ");
}  /* nodelabel */


void jointree()
{
  /* calculate the tree */
  long nc, nextnode, mini=0, minj=0, i, j, ia, ja, ii, jj, nude, iter;
  double fotu2, total, tmin, dio, djo, bi, bj, bk, dmin=0, da;
  long el[3];
  vector av;
  intvector oc;

  /* for Y. Ina speedups */
  double *R;
  R = (double *)calloc(spp,sizeof(double));

  for (i = 0; i <= spp - 2; i++) {
    for (j = i + 1; j < spp; j++) {
      da = (x[i][j] + x[j][i]) / 2.0;
      x[i][j] = da;
      x[j][i] = da;
    }
  }
  /* First initialization */
  fotu2 = spp - 2.0;
  nextnode = spp + 1;
  av = (vector)Malloc(spp*sizeof(double));
  oc = (intvector)Malloc(spp*sizeof(long));
  for (i = 0; i < spp; i++) {
    av[i] = 0.0;
    oc[i] = 1;
  }
  /* Enter the main cycle */
  if (njoin)
    iter = spp - 3;
  else
    iter = spp - 1;
  for (nc = 1; nc <= iter; nc++) {
    for (j = 2; j <= spp; j++) {
      for (i = 0; i <= j - 2; i++)
        x[j - 1][i] = x[i][j - 1];
    }
    tmin = 99999.0;
    /* Compute sij and minimize */
    if (njoin) { /* Y. Ina */
      for (i = 0; i < spp; i++) { /* Y. Ina */
          R[i] = 0.0; /* Y. Ina */
      } /* Y. Ina */
      for (ja = 2; ja <= spp; ja++) { /* Y. Ina */
        jj = enterorder[ja - 1]; /* Y. Ina */
        if (cluster[jj - 1] != NULL) { /* Y. Ina */
          for (ia = 0; ia <= ja - 2; ia++) { /* Y. Ina */
            ii = enterorder[ia]; /* Y. Ina */
            if (cluster[ii - 1] != NULL) { /* Y. Ina */
              R[ii - 1] += x[ii - 1][jj - 1]; /* Y. Ina */
              R[jj - 1] += x[ii - 1][jj - 1]; /* Y. Ina */
            } /* Y. Ina */
          } /* Y. Ina */
        } /* Y. Ina */
      } /* Y. Ina */
    } /* Y. Ina */
    for (ja = 2; ja <= spp; ja++) {
      jj = enterorder[ja - 1];
      if (cluster[jj - 1] != NULL) {
        for (ia = 0; ia <= ja - 2; ia++) {
          ii = enterorder[ia];
          if (cluster[ii - 1] != NULL) {
            if (njoin) {
      /*      diq = 0.0;
              djq = 0.0;
              dij = x[ii - 1][jj - 1];
              for (i = 0; i < spp; i++) {
                diq += x[i][ii - 1];
                djq += x[i][jj - 1];
              }
              total = fotu2 * dij - diq - djq; */ /* deleted by Y. Ina */
              total = fotu2 * x[ii - 1][jj - 1] - R[ii - 1] - R[jj - 1]; /* Y. Ina */
            } else
              total = x[ii - 1][jj - 1];
            if (total < tmin) {
              tmin = total;
              mini = ii;
              minj = jj;
            }
          }
        }
      }
    }
    /* compute lengths and print */
    if (njoin) {
      dio = 0.0;
      djo = 0.0;
      for (i = 0; i < spp; i++) {
        dio += x[i][mini - 1];
        djo += x[i][minj - 1];
      }
      dmin = x[mini - 1][minj - 1];
      dio = (dio - dmin) / fotu2;
      djo = (djo - dmin) / fotu2;
      bi = (dmin + dio - djo) * 0.5;
      bj = dmin - bi;
      bi -= av[mini - 1];
      bj -= av[minj - 1];
    } else {
      bi = x[mini - 1][minj - 1] / 2.0 - av[mini - 1];
      bj = x[mini - 1][minj - 1] / 2.0 - av[minj - 1];
      av[mini - 1] += bi;
    }
    if (progress) {
      printf("Cycle %3ld: ", iter - nc + 1);
      if (njoin)
        nodelabel(av[mini - 1] > 0.0);
      else
        nodelabel(oc[mini - 1] > 1.0);
      printf("%3ld (%10.5f) joins ", mini, bi);
      if (njoin)
        nodelabel(av[minj - 1] > 0.0);
      else
        nodelabel(oc[minj - 1] > 1.0);
      printf("%3ld (%10.5f)\n", minj, bj);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    hookup(curtree.nodep[nextnode - 1]->next, cluster[mini - 1]);
    hookup(curtree.nodep[nextnode - 1]->next->next, cluster[minj - 1]);
    cluster[mini - 1]->v = bi;
    cluster[minj - 1]->v = bj;
    cluster[mini - 1]->back->v = bi;
    cluster[minj - 1]->back->v = bj;
    cluster[mini - 1] = curtree.nodep[nextnode - 1];
    cluster[minj - 1] = NULL;
    nextnode++;
    if (njoin)
      av[mini - 1] = dmin * 0.5;
    /* re-initialization */
    fotu2 -= 1.0;
    for (j = 0; j < spp; j++) {
      if (cluster[j] != NULL) {
        if (njoin) {
          da = (x[mini - 1][j] + x[minj - 1][j]) * 0.5;
          if (mini - j - 1 < 0)
            x[mini - 1][j] = da;
          if (mini - j - 1 > 0)
            x[j][mini - 1] = da;
        } else {
          da = x[mini - 1][j] * oc[mini - 1] + x[minj - 1][j] * oc[minj - 1];
          da /= oc[mini - 1] + oc[minj - 1];
          x[mini - 1][j] = da;
          x[j][mini - 1] = da;
        }
      }
    }
    for (j = 0; j < spp; j++) {
      x[minj - 1][j] = 0.0;
      x[j][minj - 1] = 0.0;
    }
    oc[mini - 1] += oc[minj - 1];
  }
  /* the last cycle */
  nude = 1;
  for (i = 1; i <= spp; i++) {
    if (cluster[i - 1] != NULL) {
      el[nude - 1] = i;
      nude++;
    }
  }
  if (!njoin) {
    curtree.start = cluster[el[0] - 1];
    curtree.start->back = NULL;
    free(av);
    free(oc);
    return;
  }
  bi = (x[el[0] - 1][el[1] - 1] + x[el[0] - 1][el[2] - 1] - x[el[1] - 1]
        [el[2] - 1]) * 0.5;
  bj = x[el[0] - 1][el[1] - 1] - bi;
  bk = x[el[0] - 1][el[2] - 1] - bi;
  bi -= av[el[0] - 1];
  bj -= av[el[1] - 1];
  bk -= av[el[2] - 1];
  if (progress) {
    printf("last cycle:\n");
    putchar(' ');
    nodelabel(av[el[0] - 1] > 0.0);
    printf("%3ld  (%10.5f) joins ", el[0], bi);
    nodelabel(av[el[1] - 1] > 0.0);
    printf("%3ld  (%10.5f) joins ", el[1], bj);
    nodelabel(av[el[2] - 1] > 0.0);
    printf("%3ld  (%10.5f)\n", el[2], bk);
#ifdef WIN32
    phyFillScreenColor();
#endif
  }
  hookup(curtree.nodep[nextnode - 1], cluster[el[0] - 1]);
  hookup(curtree.nodep[nextnode - 1]->next, cluster[el[1] - 1]);
  hookup(curtree.nodep[nextnode - 1]->next->next, cluster[el[2] - 1]);
  cluster[el[0] - 1]->v = bi;
  cluster[el[1] - 1]->v = bj;
  cluster[el[2] - 1]->v = bk;
  cluster[el[0] - 1]->back->v = bi;
  cluster[el[1] - 1]->back->v = bj;
  cluster[el[2] - 1]->back->v = bk;
  curtree.start = cluster[el[0] - 1]->back;
  free(av);
  free(oc);
}  /* jointree */


void maketree()
{
  /* construct the tree */
  long i ;

  inputdata(replicates, printdata, lower, upper, x, reps);
  if (progress)
    putchar('\n');
  if (ith == 1)
    setuptree(&curtree, nonodes2 + 1);
  for (i = 1; i <= spp; i++)
    enterorder[i - 1] = i;
  if (jumble)
    randumize(seed, enterorder);
  for (i = 0; i < spp; i++)
    cluster[i] = curtree.nodep[i];
  jointree();
  if (outgropt)
    curtree.start = curtree.nodep[outgrno - 1]->back;
  printree(&curtree, curtree.start->next,
            treeprint, njoin, !njoin);
  if (treeprint)
    summarize();
  if (trout) {
    col = 0;
    if (njoin)
      treeout(curtree.start, &col, 0.43429448222, njoin, curtree.start);
    else
      curtree.root = curtree.start,
      treeoutr(curtree.start,&col,&curtree);
  }
  if (progress) {
    printf("\nOutput written on output file\n\n");
    if (trout)
      printf("Tree written on tree file\n\n");
  }
}  /* maketree */


int main(int argc, Char *argv[])
{  /* main program */
#ifdef MAC
  argc = 1;		/* macsetup("Neighbor","");		*/
  argv[0] = "Neighbor";
#endif
#ifdef WIN32
  phySetConsoleAttributes();
  phyClearScreen();
#endif
  openfile(&infile,INFILE,"input file", "r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file", "w",argv[0],outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  datasets = 1;
  doinit();
  if (trout)
    openfile(&outtree,OUTTREE,"output tree file", "w",argv[0],outtreename);
  ith = 1;
  while (ith <= datasets) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n",ith);
      if (progress)
        printf("Data set # %ld:\n",ith);
    }
    getinput();
    maketree();
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    ith++;
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
}





