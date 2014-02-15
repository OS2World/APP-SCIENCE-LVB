
#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1993-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxsp           4   /* maximum number of species -- must be 4 */

typedef enum {
  xx, yy, zz, ww
} simbol;

#ifndef OLDC
/* function prototypes */
void getoptions(void);
void allocrest(void);
void doinit(void);
void dnainvar_inputweights(void);
void inputoptions(void);
void dnainvar_sitecombine(void);
void makeweights(void);
void doinput(void);
void prntpatterns(void);
void makesymmetries(void);

void prntsymbol(simbol);
void prntsymmetries(void);
void tabulate(long,long,long,long,double *,double *,double *,double *);
void dnainvar_writename(long);
void writetree(long, long, long, long);
void exacttest(long, long);
void invariants(void);
void makeinv(void);
/* function prototypes */
#endif


extern FILE *infile, *outfile;
extern long spp, nonodes, endsite;
extern boolean interleaved, printdata;
extern sequence y;
extern steptr weight, alias;
extern naym *nayme;

Char infilename[100], outfilename[100];
long sites, datasets, ith;
boolean weights, progress, prntpat, printinv, mulsets, firstset;
steptr aliasweight;

long f[(long)ww - (long)xx + 1][(long)ww - (long)xx + 1]
       [(long)ww - (long)xx + 1]; /* made global from being local to makeinv */

void getoptions()
{
  /* interactively set options */
  boolean done;
  Char ch;

  fprintf(outfile, "\nNucleic acid sequence Invariants ");
  fprintf(outfile, "method, version %s\n\n",VERSION);
  putchar('\n');
  printdata = false;
  dotdiff = true;
  progress = true;
  prntpat = true;
  printinv = true;
  interleaved = true;
  do {
    if (ansi)
      printf("\033[2J\033[H");
    else
      putchar('\n');
    printf("\nNucleic acid sequence Invariants ");
    printf("method, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
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
    if (printdata)
      printf("  .  Use dot-differencing to display them  %s\n",
           dotdiff ? "Yes" : "No");
    printf("  2  Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "No"));
    printf("  3      Print out the counts of patterns");
    if (prntpat)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  4              Print out the invariants");
    if (printinv)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("\n  Y to accept these or type the letter for one to change\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if (strchr("MI01.234",ch) != NULL) {
	switch (ch) {

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

	case '.':
	  dotdiff = !dotdiff;
	  break;

        case '2':
	  progress = !progress;
    	  break;
	
	case '3':
	  prntpat = !prntpat;
	  break;

	case '4':
	  printinv = !printinv;
	  break;
	}
      } else
	printf("Not a possible option!\n");
    }
  } while (!done);
}  /* getoptions */


void allocrest()
{
  long i;

  y       = (Char **)Malloc(spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    y[i] = (Char *)Malloc(sites*sizeof(Char));
  nayme        = (naym *)Malloc(maxsp * sizeof(naym));
  weight       = (steptr)Malloc(sites * sizeof(long));
  alias        = (steptr)Malloc(sites * sizeof(long));
  aliasweight  = (steptr)Malloc(sites * sizeof(long));
}


void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &sites, &nonodes, 1);
  if (spp > maxsp){
    printf("TOO MANY SPECIES: only 4 allowed\n");
    exxit(-1);}
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, sites);
  allocrest();
}  /* doinit*/


void dnainvar_inputweights()
{
  /* input the character weights, which must be 0 or 1 */
  Char ch;
  long i;

  for (i = 1; i < nmlngth; i++)
    ch = getc(infile);
  for (i = 0; i < sites; i++) {
    do {
      if (eoln(infile)) {
	fscanf(infile, "%*[^\n]");
	getc(infile);
      }
      ch = getc(infile);
    } while (ch == ' ');
    weight[i] = 1;
    if (ch == '0')
      weight[i] = 0;
    else if (ch != '1') {
      printf("BAD WEIGHT CHARACTER: %c\n", ch);
      exxit(-1);
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  weights = true;
}  /* dnainvar_inputweights */


void inputoptions()
{
  /* input the information on the options */
  Char ch;
  long extranum, i;

  if (!firstset)
    samenumsp(&sites, ith);
  weights = false;
  extranum = 0;
  readoptions(&extranum, "W");
  for (i = 0; i < sites; i++)
    weight[i] = 1;
  for (i = 1; i <= extranum; i++) {
    matchoptions(&ch, "W");
    if (ch == 'W')
      dnainvar_inputweights();
  }
  if (weights)
    printweights(outfile, 0, sites, weight, "Sites");
}  /* inputoptions */


void dnainvar_sitecombine()
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
	tied = (tied &&
	    y[k - 1][alias[i - 1] - 1] == y[k - 1][alias[j - 1] - 1]);
	k++;
      }
      if (tied && aliasweight[j - 1] > 0) {
	aliasweight[i - 1] += aliasweight[j - 1];
	aliasweight[j - 1] = 0;
      }
      j++;
    }
    i = j - 1;
  }
}  /* dnainvar_sitecombine */


void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= sites; i++) {
    alias[i - 1] = i;
    aliasweight[i - 1] = weight[i - 1];
  }
  sitesort(sites, aliasweight);
  dnainvar_sitecombine();
  sitescrunch2(sites, 1, 2, aliasweight);
  for (i = 1; i <= sites; i++) {
    weight[i - 1] = aliasweight[i - 1];
    if (weight[i - 1] > 0)
      endsite = i;
  }
}  /* makeweights */


void doinput()
{
  /* reads the input data */
  inputoptions();
  inputdata(sites);
  makeweights();
}  /* doinput */



void prntpatterns()
{
  /* print out patterns */
  long i, j;

  fprintf(outfile, "\n   Pattern");
  if (prntpat)
    fprintf(outfile, "   Number of times");
  fprintf(outfile, "\n\n");
  for (i = 0; i < endsite; i++) {
    fprintf(outfile, "     ");
    for (j = 0; j < spp; j++)
      putc(y[j][alias[i] - 1], outfile);
    if (prntpat)
      fprintf(outfile, "  %8ld", weight[i]);
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* prntpatterns */


void makesymmetries()
{
  /* get frequencies of symmetrized patterns */
  long i, j;
  boolean drop, usedz;
  Char ch, ch1, zchar;
  simbol s1, s2, s3;
  simbol t[maxsp - 1];

  for (s1 = xx; (long)s1 <= (long)ww; s1 = (simbol)((long)s1 + 1)) {
    for (s2 = xx; (long)s2 <= (long)ww; s2 = (simbol)((long)s2 + 1)) {
      for (s3 = xx; (long)s3 <= (long)ww; s3 = (simbol)((long)s3 + 1))
	f[(long)s1 - (long)xx][(long)s2 - (long)xx]
	  [(long)s3 - (long)xx] = 0;
    }
  }
  for (i = 0; i < endsite; i++) {
    drop = false;
    for (j = 0; j < spp; j++) {
      ch = y[j][alias[i] - 1];
      drop = (drop ||
	      (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T' && ch != 'U'));
    }
    ch1 = y[0][alias[i] - 1];
    if (!drop) {
      usedz = false;
      zchar = ' ';
      for (j = 2; j <= spp; j++) {
	ch = y[j - 1][alias[i] - 1];
	if (ch == ch1)
	  t[j - 2] = xx;
	else if ((ch1 == 'A' && ch == 'G') || (ch1 == 'G' && ch == 'A') ||
		 (ch1 == 'C' && (ch == 'T' || ch == 'U')) ||
		 ((ch1 == 'T' || ch1 == 'U') && ch == 'C'))
	  t[j - 2] = yy;
	else if (!usedz) {
	  t[j - 2] = zz;
	  usedz = true;
	  zchar = ch;
	} else if (usedz && ch == zchar)
	  t[j - 2] = zz;
	else if (usedz && ch != zchar)
	  t[j - 2] = ww;
      }
      f[(long)t[0] - (long)xx][(long)t[1] - (long)xx]
	[(long)t[2] - (long)xx] += weight[i];
    }
  }
}  /* makesymmetries */


void prntsymbol(simbol s)
{
  /* print 1, 2, 3, 4 as appropriate */
  switch (s) {

  case xx:
    putc('1', outfile);
    break;

  case yy:
    putc('2', outfile);
    break;

  case zz:
    putc('3', outfile);
    break;

  case ww:
    putc('4', outfile);
    break;
  }
}  /* prntsymbol */


void prntsymmetries()
{
  /* print out symmetrized pattern numbers */
  simbol s1, s2, s3;

  fprintf(outfile, "\nSymmetrized patterns (1, 2 = the two purines  ");
  fprintf(outfile, "and  3, 4 = the two pyrimidines\n");
  fprintf(outfile, "                  or  1, 2 = the two pyrimidines  ");
  fprintf(outfile, "and  3, 4 = the two purines)\n\n");
  for (s1 = xx; (long)s1 <= (long)ww; s1 = (simbol)((long)s1 + 1)) {
    for (s2 = xx; (long)s2 <= (long)ww; s2 = (simbol)((long)s2 + 1)) {
      for (s3 = xx; (long)s3 <= (long)ww; s3 = (simbol)((long)s3 + 1)) {
	if (f[(long)s1 - (long)xx][(long)s2 - (long)xx]
	    [(long)s3 - (long)xx] > 0) {
	  fprintf(outfile, "     1");
	  prntsymbol(s1);
	  prntsymbol(s2);
	  prntsymbol(s3);
	  if (prntpat)
	    fprintf(outfile, "   %7ld",
		    f[(long)s1 - (long)xx][(long)s2 - (long)xx]
		    [(long)s3 - (long)xx]);
	  putc('\n', outfile);
	}
      }
    }
  }
}  /* prntsymmetries */


void tabulate(long mm, long nn, long pp, long qq, double *mr,
			double *nr, double *pr, double *qr)
{
  /* make quadratic invariant, table, chi-square */
  long total;
  double k, TEMP;

  fprintf(outfile, "\n   Contingency Table\n\n");
  fprintf(outfile, "%7ld%6ld\n", mm, nn);
  fprintf(outfile, "%7ld%6ld\n\n", pp, qq);
  *mr = (long)(mm);
  *nr = (long)(nn);
  *pr = (long)pp;
  *qr = (long)qq;
  total = mm + nn + pp + qq;
  if (printinv)
    fprintf(outfile, "   Quadratic invariant = %15.1f\n\n",
	    (*nr) * (*pr) - (*mr) * (*qr));
  fprintf(outfile, "   Chi-square = ");
  TEMP = (*mr) * (*qr) - (*nr) * (*pr);
  k = total * (TEMP * TEMP) / (((*mr) + (*nr)) * ((*mr) + (*pr)) *
			       ((*nr) + (*qr)) * ((*pr) + (*qr)));
  fprintf(outfile, "%10.5f", k);
  if ((*mr) * (*qr) > (*nr) * (*pr) && k > 2.71)
    fprintf(outfile, " (P < 0.05)\n");
  else
    fprintf(outfile, " (not significant)\n");
  fprintf(outfile, "\n\n");
}  /* tabulate */


void dnainvar_writename(long m)
{
  /* write out a species name */
  long i, n;

  n = nmlngth;
  while (nayme[m - 1][n - 1] == ' ')
    n--;
  if (n == 0)
    n = 1;
  for (i = 0; i < n; i++)
    putc(nayme[m - 1][i], outfile);
}  /* dnainvar_writename */


void writetree(long i, long j, long k, long l)
{
  /* write out tree topology ((i,j),(k,l)) using names */
  fprintf(outfile, "((");
  dnainvar_writename(i);
  putc(',', outfile);
  dnainvar_writename(j);
  fprintf(outfile, "),(");
  dnainvar_writename(k);
  putc(',', outfile);
  dnainvar_writename(l);
  fprintf(outfile, "))\n");
}  /* writetree */


void exacttest(long m, long n)
{
  /* exact binomial test that m <= n */
  long i;
  double p, sum;

  p = 1.0;
  for (i = 1; i <= m + n; i++)
    p /= 2.0;
  sum = p;
  for (i = 1; i <= n; i++) {
    p = p * (m + n - i + 1) / i;
    sum += p;
  }
  fprintf(outfile, "      %7.4f", sum);
  if (sum <= 0.05)
    fprintf(outfile, "              yes\n");
  else
    fprintf(outfile, "               no\n");
}  /* exacttest */


void invariants()
{
  /* compute invariants */
  long  m, n, p, q;
  double L1, L2, L3;
  double mr,nr,pr,qr;

  fprintf(outfile, "\nTree topologies (unrooted): \n\n");
  fprintf(outfile, "    I:  ");
  writetree(1, 2, 3, 4);
  fprintf(outfile, "   II:  ");
  writetree(1, 3, 2, 4);
  fprintf(outfile, "  III:  ");
  writetree(1, 4, 2, 3);
  fprintf(outfile, "\n\nLake's linear invariants\n");
  fprintf(outfile,
    " (these are expected to be zero for the two incorrect tree topologies.\n");
  fprintf(outfile,
	  "  This is tested by testing the equality of the two parts\n");
  fprintf(outfile,
	  "  of each expression using a one-sided exact binomial test.\n");
  fprintf(outfile,
    "  The null hypothesis is that the first part is no larger than the second.)\n\n");
  fprintf(outfile, " Tree                           ");
  fprintf(outfile, "  Exact test P value    Significant?\n\n");
  m = f[(long)yy - (long)xx][(long)zz - (long)xx]
      [(long)ww - (long)xx] + f[0][(long)zz - (long)xx]
      [(long)zz - (long)xx];
  n = f[(long)yy - (long)xx][(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[0][(long)zz - (long)xx]
      [(long)ww - (long)xx];
  fprintf(outfile, "   I  %5ld    - %5ld   = %5ld", m, n, m - n);
  exacttest(m, n);
  m = f[(long)zz - (long)xx][(long)yy - (long)xx]
      [(long)ww - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)zz - (long)xx];
  n = f[(long)zz - (long)xx][(long)yy - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)ww - (long)xx];
  fprintf(outfile, "   II %5ld    - %5ld   = %5ld", m, n, m - n);
  exacttest(m, n);
  m = f[(long)zz - (long)xx][(long)ww - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx][0];
  n = f[(long)zz - (long)xx][(long)zz - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)ww - (long)xx][0];
  fprintf(outfile, "   III%5ld    - %5ld   = %5ld", m, n, m - n);
  exacttest(m, n);
  fprintf(outfile, "\n\nCavender's quadratic invariants (type L)");
  fprintf(outfile, " using purines vs. pyrimidines\n");
  fprintf(outfile,
	  " (these are expected to be zero, and thus have a nonsignificant\n");
  fprintf(outfile, "  chi-square, for the correct tree topology)\n");
  fprintf(outfile, "They will be misled if there are substantially\n");
  fprintf(outfile, "different evolutionary rate between sites, or\n");
  fprintf(outfile, "different purine:pyrimidine ratios from 1:1.\n\n");
  fprintf(outfile, "  Tree I:\n");
  m = f[0][0][0] + f[0][(long)yy - (long)xx]
      [(long)yy - (long)xx] + f[0][(long)zz - (long)xx]
      [(long)zz - (long)xx];
  n = f[0][0][(long)yy - (long)xx] + f[0][0]
      [(long)zz - (long)xx] + f[0][(long)yy - (long)xx][0] + f[0]
      [(long)yy - (long)xx][(long)zz - (long)xx] + f[0]
      [(long)zz - (long)xx][0] + f[0][(long)zz - (long)xx]
      [(long)yy - (long)xx] + f[0][(long)zz - (long)xx]
      [(long)ww - (long)xx];
  p = f[(long)yy - (long)xx][0][0] + f[(long)yy - (long)xx]
      [(long)yy - (long)xx]
      [(long)yy - (long)xx] + f[(long)yy - (long)xx]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx][0]
      [0] + f[(long)zz - (long)xx][(long)yy - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)ww - (long)xx][(long)ww - (long)xx];
  q = f[(long)yy - (long)xx][0][(long)yy - (long)xx] +
      f[(long)yy - (long)xx][0][(long)zz - (long)xx] +
      f[(long)yy - (long)xx][(long)yy - (long)xx][0] +
      f[(long)yy - (long)xx][(long)yy - (long)xx][(long)zz - (long)xx] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][0] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][(long)yy - (long)xx] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][0][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][0][(long)zz - (long)xx] +
      f[(long)zz - (long)xx][0][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][0] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][(long)zz - (long)xx] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][0] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][0] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][(long)zz - (long)xx];

  nr = n;
  pr = p;
  mr = m;
  qr = q;
  L1 = nr * pr - mr * qr;
  tabulate(m, n, p, q, &mr,&nr,&pr,&qr);
  fprintf(outfile, "  Tree II:\n");
  m = f[0][0][0] + f[(long)yy - (long)xx][0]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)zz - (long)xx];
  n = f[0][0][(long)yy - (long)xx] + f[0][0]
      [(long)zz - (long)xx] + f[(long)yy - (long)xx][0]
      [0] + f[(long)yy - (long)xx][0]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx][0]
      [0] + f[(long)zz - (long)xx][0]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)ww - (long)xx];
  p = f[0][(long)yy - (long)xx][0] + f[(long)yy - (long)xx]
      [(long)yy - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)yy - (long)xx][(long)zz - (long)xx] + f[0]
      [(long)zz - (long)xx][0] + f[(long)yy - (long)xx]
      [(long)zz - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)ww - (long)xx][(long)zz - (long)xx];
  q = f[0][(long)yy - (long)xx][(long)yy - (long)xx] + f[0]
      [(long)yy - (long)xx][(long)zz - (long)xx] +
      f[(long)yy - (long)xx][(long)yy - (long)xx][0] +
      f[(long)yy - (long)xx][(long)yy - (long)xx][(long)zz - (long)xx] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][0] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][(long)ww - (long)xx] +
      f[0][(long)zz - (long)xx][(long)yy - (long)xx] + f[0]
      [(long)zz - (long)xx][(long)zz - (long)xx] + f[0]
      [(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][0] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][(long)zz - (long)xx] +
      f[(long)yy - (long)xx][(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][0] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][0] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][(long)ww - (long)xx];
  nr = n;
  pr = p;
  mr = m;
  qr = q;
  L2 = nr * pr - mr * qr;
  tabulate(m, n, p, q, &mr,&nr,&pr,&qr);
  fprintf(outfile, "  Tree III:\n");
  m = f[0][0][0] + f[(long)yy - (long)xx][(long)yy - (long)xx]
      [0] + f[(long)zz - (long)xx][(long)zz - (long)xx][0];
  n = f[(long)yy - (long)xx][0][0] + f[(long)zz - (long)xx][0]
      [0] + f[0][(long)yy - (long)xx][0] + f[(long)zz - (long)xx]
      [(long)yy - (long)xx][0] + f[0][(long)zz - (long)xx]
      [0] + f[(long)yy - (long)xx][(long)zz - (long)xx]
      [0] + f[(long)zz - (long)xx][(long)ww - (long)xx][0];
  p = f[0][0][(long)yy - (long)xx] + f[(long)yy - (long)xx]
      [(long)yy - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx][(long)yy - (long)xx] + f[0][0]
      [(long)zz - (long)xx] + f[(long)yy - (long)xx]
      [(long)yy - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)zz - (long)xx][(long)ww - (long)xx];
  q = f[(long)yy - (long)xx][0][(long)yy - (long)xx] + f[(long)zz - (long)xx]
      [0][(long)yy - (long)xx] + f[0][(long)yy - (long)xx][(long)yy - (long)xx] +
      f[(long)zz - (long)xx][(long)yy - (long)xx][(long)yy - (long)xx] +
      f[0][(long)zz - (long)xx]
      [(long)yy - (long)xx] + f[(long)yy - (long)xx][(long)zz - (long)xx]
      [(long)yy - (long)xx] + f[(long)zz - (long)xx][(long)ww - (long)xx]
      [(long)yy - (long)xx] + f[(long)yy - (long)xx][0]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)zz - (long)xx] + f[0][(long)zz - (long)xx]
      [(long)ww - (long)xx] + f[0][(long)yy - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)yy - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)yy - (long)xx][(long)ww - (long)xx] + f[0]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)yy - (long)xx]
      [(long)zz - (long)xx]
      [(long)zz - (long)xx] + f[(long)zz - (long)xx]
      [(long)ww - (long)xx]
      [(long)ww - (long)xx] + f[(long)zz - (long)xx][0]
      [(long)ww - (long)xx] + f[(long)yy - (long)xx]
      [(long)zz - (long)xx][(long)ww - (long)xx] +
      f[(long)zz - (long)xx][(long)ww - (long)xx][(long)zz - (long)xx];
  nr = n;
  pr = p;
  mr = m;
  qr = q;
  L3 = nr * pr - mr * qr;
  tabulate(m, n, p, q, &mr,&nr,&pr,&qr);
  fprintf(outfile, "\n\nCavender's quadratic invariants (type K)");
  fprintf(outfile, " using purines vs. pyrimidines\n");
  fprintf(outfile,
	  " (these are expected to be zero for the correct tree topology)\n");
  fprintf(outfile, "They will be misled if there are substantially\n");
  fprintf(outfile, "different evolutionary rate between sites, or\n");
  fprintf(outfile, "different purine:pyrimidine ratios from 1:1.\n");
  fprintf(outfile, "No statistical test is done on them here.\n\n");
  fprintf(outfile, "  Tree I:   %15.1f\n", L2 - L3);
  fprintf(outfile, "  Tree II:  %15.1f\n", L3 - L1);
  fprintf(outfile, "  Tree III: %15.1f\n\n", L1 - L2);
}  /* invariants */


void makeinv()
{
  /* print out patterns and compute invariants */

  prntpatterns();
  makesymmetries();
  prntsymmetries();
  if (printinv)
    invariants();
}  /* makeinv */


int main(int argc, Char *argv[])
{  /* DNA Invariants */
#ifdef MAC
  argc = 1;		/* macsetup("Dnainvar","");		*/
  argv[0] = "Dnainvar";         
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
  firstset = true;
  datasets = 1;
  doinit();
  for (ith = 1; ith <= datasets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      if (progress)
        printf("\nData set # %ld:\n",ith);
      fprintf(outfile, "Data set # %ld:\n\n",ith);
    }
    makeinv();
  }
  if (progress) {
    putchar('\n');
    printf("Output written to output file\n");
    putchar('\n');
  }
  FClose(outfile);
  FClose(infile);
#ifdef MAC
  fixmacfile(outfilename);
#endif
  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* DNA Invariants */
