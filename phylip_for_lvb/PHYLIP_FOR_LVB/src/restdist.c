
#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1994-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define initialv        0.1     /* starting value of branch length          */
#define iterationsr     20      /* how many Newton iterations per distance  */

extern FILE *infile, *outfile;
extern long spp, endsite;
extern boolean interleaved, printdata;
extern sequence y;
extern steptr weight, alias;
extern naym *nayme;

#ifndef OLDC
/* function prototypes */
void restdist_inputnumbers(void);
void getoptions(void);
void allocrest(void);
void doinit(void);
void inputoptions(void);
void restdist_inputdata(void);
void restdist_sitesort(void);
void restdist_sitecombine(void);
void makeweights(void);
void makev(long, long, double *);
void makedists(void);
void writedists(void);
void getinput(void);
/* function prototypes */
#endif


Char infilename[100], outfilename[100]; 
long sites, weightsum, datasets, ith;
boolean  restsites, neili, gama, weights, lower,
           progress, mulsets, firstset;
double ttratio, fracchange, cvi, sitelength, xi, xv;
double **d;
steptr aliasweight;
Char progname[20];
Char ch;

void restdist_inputnumbers()
{
  /* read and print out numbers of species and sites */
  fscanf(infile, "%ld%ld", &spp, &sites);
}  /* restdist_inputnumbers */

void getoptions()
{
  /* interactively set options */
  Char ch;

  putchar('\n');
  sitelength = 6.0;
  neili = false;
  gama = false;
  cvi = 0.0;
  weights = false;
  lower = false;
  printdata = false;
  progress = true;
  restsites = true;
  interleaved = true;
  ttratio = 2.0;
  for (;;) {
    printf("%s", (ansi ? ("\033[2J\033[H") : "\n"));
    printf("\nRestriction site or fragment distances, ");
    printf("version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  R           Restriction sites or fragments?  %s\n",
           (restsites ? "Sites" : "Fragments"));
    if (!restsites)
      printf("  N        Original or modified Nei/Li model?  %s\n",
             (neili ? "Original" : "Modified"));
    if (restsites || !neili) {
      printf("  G  Gamma distribution of rates among sites?");
      if (!gama)
        printf("  No\n");
      else
        printf("  Yes\n");
      printf("  T            Transition/transversion ratio?  %f\n", ttratio);
    }
    printf("  S                              Site length? %4.1f\n",sitelength);
    printf("  L                  Form of distance matrix?  %s\n",
           (lower ? "Lower-triangular" : "Square"));
    printf("  M               Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  I              Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0       Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1       Print out the data at start of run?  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2     Print indications of progress of run?  %s\n",
           (progress ? "Yes" : "No"));
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
    if (strchr("RDNGTSLMI012",ch) != NULL){
      switch (ch) {
        
      case 'R':
        restsites = !restsites;
        break;
        
      case 'G':
        if (restsites || !neili)
          gama = !gama;
        break;

      case 'N':
        if (!restsites)
          neili = !neili;
        break;

      case 'T':
        if (restsites || !neili)
          initratio(&ttratio);
        break;

      case 'S':
        do {
          printf("New Sitelength?\n");
#ifdef WIN32
          phyFillScreenColor();
#endif
          scanf("%lf%*[^\n]", &sitelength);
          getchar();
          if (sitelength < 1.0)
            printf("BAD RESTRICTION SITE LENGTH: %f\n", sitelength); 
        } while (sitelength < 1.0);
        break;
        
      case 'L':
        lower = !lower;
        break;

      case 'M':
        mulsets = !mulsets;
        if (mulsets)
          initdatasets(&datasets);
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
        
      }
    } else
      printf("Not a possible option!\n");
  }
  if (gama) {
    do {
      printf(
"\nCoefficient of variation of substitution rate among sites (must be positive)?\n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      scanf("%lf%*[^\n]", &cvi);
      getchar();
    } while (cvi <= 0.0);
    cvi = 1.0 / (cvi * cvi);
    printf("\n");
  }
  xi = (ttratio - 0.5)/(ttratio + 0.5);
  xv = 1.0 - xi;
  fracchange = xi*0.5 + xv*0.75;
}  /* getoptions */

void allocrest()
{
  long i;

  y = (Char **)Malloc(spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    y[i] = (Char *)Malloc(sites*sizeof(Char));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  weight = (steptr)Malloc((sites+1)*sizeof(long));
  alias = (steptr)Malloc((sites+1)*sizeof(long));
  aliasweight = (steptr)Malloc((sites+1)*sizeof(long));
  d = (double **)Malloc(spp*sizeof(double *));
  for (i = 0; i < spp; i++)
    d[i] = (double*)Malloc(spp*sizeof(double));
}  /* allocrest */


void doinit()
{
  /* initializes variables */
  restdist_inputnumbers();
  getoptions();
  if (printdata)
    fprintf(outfile, "\n %4ld Species, %4ld Sites\n", spp, sites);
  allocrest();
}  /* doinit */


void inputoptions()
{
  /* read the options information */
  Char ch;
  long i, extranum, cursp, curst;

  if (!firstset) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    fscanf(infile, "%ld%ld", &cursp, &curst);
    if (cursp != spp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n",
             ith);
      exxit(-1);
    }
    sites = curst;
  }
  for (i = 1; i <= sites; i++)
    weight[i] = 1;
  weightsum = sites;
  extranum = 0;
  fscanf(infile, "%*[ 0-9]");
  readoptions(&extranum, "W");
  for (i = 1; i <= extranum; i++) {
      matchoptions(&ch, "W");
      inputweights2(1, sites+1, &weightsum, weight, &weights, "RESTDIST");
  }
}  /* inputoptions */


void restdist_inputdata()
{
  /* read the species and sites data */
  long i, j, k, l, sitesread = 0, sitesnew = 0;
  Char ch;
  boolean allread, done;

  if (printdata)
    putc('\n', outfile);
  j = nmlngth + (sites + (sites - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 39)
    j = 39;
  if (printdata) {
    fprintf(outfile, "Name");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "Sites\n");
    fprintf(outfile, "----");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "-----\n\n");
  }
  sitesread = 0;
  allread = false;
  while (!(allread)) {
    allread = true;
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    i = 1;
    while (i <= spp ) {
      if ((interleaved && sitesread == 0) || !interleaved)
        initname(i - 1);
      if (interleaved)
        j = sitesread;
      else
        j = 0;
      done = false;
      while (!done && !eoff(infile)) {
        if (interleaved)
          done = true;
        while (j < sites && !(eoln(infile) || eoff(infile))) {
          ch = getc(infile);
          if (ch == '\n')
            ch = ' ';
          if (ch == ' ')
            continue;
          uppercase(&ch);
          if (ch != '1' && ch != '0' && ch != '+' && ch != '-' && ch != '?') {
            printf(" WARNING -- BAD SYMBOL %c",ch);
            printf(" AT POSITION %5ld OF SPECIES %3ld\n",j,i);
            exxit(-1);
          }
          if (ch == '1')
            ch = '+';
          if (ch == '0')
            ch = '-';
          j++;
          y[i - 1][j - 1] = ch;
        }
        if (interleaved)
          continue;
        if (j < sites) {
          fscanf(infile, "%*[^\n]");
          getc(infile);
        } else if (j == sites)
          done = true;
      }
      if (interleaved && i == 1)
        sitesnew = j;
      fscanf(infile, "%*[^\n]");
      getc(infile);
      if ((interleaved && j != sitesnew ) || ((!interleaved) && j != sites)){
        printf("ERROR: SEQUENCES OUT OF ALIGNMENT\n");
        exxit(-1);}
      i++;
    }
    if (interleaved) {
      sitesread = sitesnew;
      allread = (sitesread == sites);
    } else
      allread = (i > spp);
  }
  if (printdata) {
    for (i = 1; i <= ((sites - 1) / 60 + 1); i++) {
      for (j = 0; j < spp; j++) {
        for (k = 0; k < nmlngth; k++)
          putc(nayme[j][k], outfile);
        fprintf(outfile, "   ");
        l = i * 60;
        if (l > sites)
          l = sites;
        for (k = (i - 1) * 60 + 1; k <= l; k++) {
          putc(y[j][k - 1], outfile);
          if (k % 10 == 0 && k % 60 != 0)
            putc(' ', outfile);
        }
        putc('\n', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
}  /* restdist_inputdata */


void restdist_sitesort()
{
  /* Shell sort keeping alias, aliasweight in same order */
  long gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = sites / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= sites; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip) {
        jj = alias[j];
        jg = alias[j + gap];
        flip = false;
        tied = true;
        k = 1;
        while (k <= spp && tied) {
          flip = (y[k - 1][jj - 1] > y[k - 1][jg - 1]);
          tied = (tied && y[k - 1][jj - 1] == y[k - 1][jg - 1]);
          k++;
        }
        if (tied) {
          aliasweight[j] += aliasweight[j + gap];
          aliasweight[j + gap] = 0;
        }
        if (!flip)
          break;
        itemp = alias[j];
        alias[j] = alias[j + gap];
        alias[j + gap] = itemp;
        itemp = aliasweight[j];
        aliasweight[j] = aliasweight[j + gap];
        aliasweight[j + gap] = itemp;
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* restdist_sitesort */


void restdist_sitecombine()
{
  /* combine sites that have identical patterns */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < sites) {
    j = i + 1;
    tied = true;
    while (j <= sites && tied) {
      k = 1;
      while (k <= spp && tied) {
        tied = (tied && y[k - 1][alias[i] - 1] == y[k - 1][alias[j] - 1]);
        k++;
      }
      if (tied && aliasweight[j] > 0) {
        aliasweight[i] += aliasweight[j];
        aliasweight[j] = 0;
        alias[j] = alias[i];
      }
      j++;
    }
    i = j - 1;
  }
}  /* restdist_sitecombine */


void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= sites; i++) {
    alias[i] = i;
    aliasweight[i] = weight[i];
  }
  restdist_sitesort();
  restdist_sitecombine();
  sitescrunch2(sites + 1, 2, 3, aliasweight);
  for (i = 1; i <= sites; i++) {
    weight[i] = aliasweight[i];
    if (weight[i] > 0)
      endsite = i;
  }
  weight[0] = 1;
}  /* makeweights */


void makev(long m, long n, double *v)
{
  /* compute one distance */
  long i, ii, it, numerator, denominator;
  double f, g=0, h, b, c, p1, p2, p3, q1, pp, tt, delta, vv, slope;

  numerator = 0;
  denominator = 0;
  for (i = 0; i < endsite; i++) {
    ii = alias[i];
    if ((y[m-1][ii-1] == '+') || (y[n-1][ii-1] == '+')) {
      denominator += weight[i];
      if ((y[m-1][ii-1] == '+') && (y[n-1][ii-1] == '+')) {
        numerator += weight[i];
      }
    }
  }
  f = 2*numerator/(double)(denominator+numerator);
  if (restsites) {
    if (exp(-sitelength*1.38629436) > f) {
      printf("\nWARNING: INFINITE DISTANCE BETWEEN ");
      printf(" SPECIES %3ld AND %3ld\n", m, n);
      exxit(-1);
    }
  }
  if (!restsites) {
    if (!neili) {
      f = (sqrt(f*(f+8.0))-f)/2.0;
    }
    else {
      g = initialv;
      h = g;
      delta = g;
      it = 0;
      while (fabs(delta) > 0.00002 && it < iterationsr) {
        it++;
        h = g;
        g = exp(0.25*log(f * (3-2*g)));
        delta = g - h;
      }
    }
  }
  if ((!restsites) && neili)
    vv = - (2.0/sitelength) * log(g);
  else {
    pp = exp(log(f)/sitelength);
    delta = initialv;
    tt = delta;
    it = 0;
    while (fabs(delta) > 0.00002 && it < iterationsr) {
      it++;
      if (gama) {
        p1 = exp(-cvi * log(1 + tt / cvi));
        p2 = exp(-cvi * log(1 + xv * tt / cvi))
              - exp(-cvi * log(1 + tt / cvi));
        p3 = 1.0 - exp(-cvi * log(1 + xv * tt / cvi));
      } else {
        p1 = exp(-tt);
        p2 = exp(-xv * tt) - exp(-tt);
        p3 = 1.0 - exp(-xv * tt);
      }
      q1 = p1 + p2 / 2.0 + p3 / 4.0;
      g = q1 - pp;
      if (gama)
        slope = - 0.5 * (1 / (1 + tt / cvi)) * exp(-cvi * log(1 + tt / cvi))
                  - 0.25 * (xv / (1 + xv * tt / cvi)) *
                  exp(-cvi * log(1 + xv * tt / cvi));
      else
        slope = - 0.5*exp(-tt) - 0.25*exp(-xv*tt);
      if (g < 0.0)
        delta = fabs(delta) / -2.0;
      else
        delta = fabs(delta);
      tt += delta;
    }
    vv = fracchange * tt;
  }
  *v = vv;
}  /* makev */


void makedists()
{
  /* compute distance matrix */
  long i, j;
  double v;

  if (progress) {
    printf("Distances calculated for species\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
  }
  for (i = 0; i < spp; i++)
    d[i][i] = 0.0;
  for (i = 1; i < spp; i++) {
    if (progress) {
      printf("    ");
      for (j = 0; j < nmlngth; j++)
        putchar(nayme[i - 1][j]);
      printf("   ");
    }
    for (j = i + 1; j <= spp; j++) {
      makev(i, j, &v);
      d[i - 1][j - 1] = v;
      d[j - 1][i - 1] = v;
      if (progress)
        putchar('.');
    }
    if (progress) {
      putchar('\n');
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
  }
  if (progress) {
    printf("    ");
    for (j = 0; j < nmlngth; j++)
      putchar(nayme[spp - 1][j]);
    putchar('\n');
#ifdef WIN32
    phyFillScreenColor();
#endif
  }
}  /* makedists */


void writedists()
{
  /* write out distances */
  long i, j, k;

  if (!printdata)
    fprintf(outfile, "%5ld\n", spp);
  for (i = 0; i < spp; i++) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[i][j], outfile);
    if (lower)
      k = i;
    else
      k = spp;
    for (j = 1; j <= k; j++) {
      fprintf(outfile, "%8.4f", d[i][j - 1]);
      if ((j + 1) % 9 == 0 && j < k)
        putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  if (progress)
    printf("\nDistances written to file\n\n");
}  /* writedists */


void getinput()
{
  /* reads the input data */
  inputoptions();
  restdist_inputdata();
  makeweights();
}  /* getinput */


int main(int argc, Char *argv[])
{  /* distances from restriction sites or fragments */
#ifdef MAC
  argc = 1;		/* macsetup("Restdist",""); */
  argv[0] = "Restdist";
#endif
#ifdef WIN32
  phySetConsoleAttributes();
  phyClearScreen();
#endif
  strcpy(progname,argv[0]);
  openfile(&infile,INFILE,"input data file","r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file","w",argv[0],outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  datasets = 1;
  firstset = true;
  doinit();
  for (ith = 1; ith <= datasets; ith++) {
    getinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1 && progress)
      printf("\nData set # %ld:\n\n",ith);
    makedists();
    writedists();
  }
  FClose(infile);
  FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* distances from restriction sites or fragments */
