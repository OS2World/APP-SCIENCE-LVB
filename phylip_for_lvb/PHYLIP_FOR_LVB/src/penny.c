
#include "phylip.h"
#include "disc.h"
#include "wagner.h"

/* version 3.6. (c) Copyright 1993-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxtrees        100  /* maximum number of trees to be printed out   */
#define often           100  /* how often to notify how many trees examined */
#define many            1000 /* how many multiples of howoften before stop  */

typedef long *treenumbers;
typedef double *valptr;
typedef long *placeptr;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   inputoptions(void);
void   doinput(void);
void   supplement(bitptr);
void   evaluate(node2 *);
void   addtraverse(node2 *,node2 *,node2 *,long *,long *,valptr,placeptr);
void   addit(long);
void   reroot(node2 *);
void   describe(void);
void   maketree(void);
/* function prototypes */
#endif


extern FILE *infile, *outfile, *outtree;
extern long spp ,chars, words, bits, nextree;
extern boolean printdata;
extern steptr extras, weight;
extern naym *nayme;

Char infilename[100], outfilename[100], outtreename[100];
node2 *root;
long outgrno, rno, howmany, howoften, col, datasets, ith;
/*  outgrno indicates outgroup					*/

extern boolean ibmpc, ansi;
boolean weights, thresh, ancvar, questions, allsokal, allwagner,
	       mixture, simple, trout, noroot, didreroot, outgropt,
	       progress, treeprint, stepbox, ancseq, mulsets, firstset;
boolean *ancone, *anczero, *ancone0, *anczero0;
pointptr2 treenode;         /* pointers to all nodes in tree   */
double fracdone, fracinc;
double threshold;
double *threshwt;
bitptr wagner, wagner0;
boolean *added;
Char *guess;
steptr numsteps;
long **bestorders, **bestrees;
steptr numsone, numszero;
gbit *garbage;

long examined, mults;
boolean firsttime, done, full;
double like, bestyet;
treenumbers current, order;
long fullset;
bitptr zeroanc, oneanc;
bitptr suppsteps;


void getoptions()
{
  /* interactively set options */
  Char ch;

  fprintf(outfile, "\nPenny algorithm, version %s\n",VERSION);
  fprintf(outfile, " branch-and-bound to find all");
  fprintf(outfile, " most parsimonious trees\n\n");
  howoften = often;
  howmany = many;
  outgrno = 1;
  outgropt = false;
  simple = true;
  thresh = false;
  threshold = spp;
  trout = true;
  weights = false;
  ancvar = false;
  allsokal = false;
  allwagner = true;
  mixture = false;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  for(;;) {
    printf(ansi ? "\033[2J\033[H" : "\n");
    printf("\nPenny algorithm, version %s\n",VERSION);
    printf(" branch-and-bound to find all most parsimonious trees\n\n");
    printf("Settings for this run:\n");
    printf("  X                     Use Mixed method?  %s\n",
	   mixture ? "Yes" : "No");
    printf("  P                     Parsimony method?  %s\n",
	   (allwagner && !mixture)  ? "Wagner"       :
	   (!(allwagner || mixture)) ? "Camin-Sokal"  : "(methods in mixture)");
    printf("  F        How often to report, in trees:%5ld\n",howoften);
    printf("  H        How many groups of%5ld trees:%6ld\n",howoften,howmany);
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at species number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
    printf("  S           Branch and bound is simple?  %s\n",
	   simple ? "Yes" : "No. reconsiders order of species");
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per char.\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  A   Use ancestral states in input file?  %s\n",
	   ancvar ? "Yes" : "No");
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
	   ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
	   printdata ? "Yes" : "No");
    printf("  2  Print indications of progress of run  %s\n",
	   progress ? "Yes" : "No");
    printf("  3                        Print out tree  %s\n",
	   treeprint ? "Yes" : "No");
    printf("  4     Print out steps in each character  %s\n",
	   stepbox ? "Yes" : "No");
    printf("  5     Print states at all nodes of tree  %s\n",
	   ancseq ? "Yes" : "No");
    printf("  6       Write out trees onto tree file?  %s\n",
	   trout ? "Yes" : "No");
    printf("\nAre these settings correct?");
    printf(" (type Y or the letter for one to change)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("HFSOMPATX1234560",ch) != NULL){
      switch (ch) {
	
      case 'X':
	mixture = !mixture;
	break;
	
      case 'P':
	allwagner = !allwagner;
	break;
	
      case 'A':
	ancvar = !ancvar;
	break;
	
      case 'H':
        inithowmany(&howmany, howoften);
	break;
	
      case 'F':
        inithowoften(&howoften);
	break;
	
      case 'S':
	simple = !simple;
	break;
	
      case 'O':
	outgropt = !outgropt;
	if (outgropt)
	  initoutgroup(&outgrno, spp);
	else outgrno = 1;
	break;
	
      case 'T':
	thresh = !thresh;
	if (thresh)
	  initthreshold(&threshold);
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
	stepbox = !stepbox;
	break;
	
      case '5':
	ancseq = !ancseq;
	break;
	
      case '6':
	trout = !trout;
	break;
      }
    } else
      printf("Not a possible option!\n");
  }
  allsokal = (!allwagner && !mixture);
}  /* getoptions */


void allocrest()
{
  long i;

  weight = (steptr)Malloc(chars*sizeof(steptr));
  threshwt = (double *)Malloc(chars*sizeof(double));
  bestorders = (long **)Malloc(maxtrees*sizeof(long *));
  bestrees = (long **)Malloc(maxtrees*sizeof(long *));
  for (i = 1; i <= maxtrees; i++) {
    bestorders[i - 1] = (long *)Malloc(spp*sizeof(long));
    bestrees[i - 1] = (long *)Malloc(spp*sizeof(long));
  }
  numsteps = (steptr)Malloc(chars*sizeof(steptr));
  guess = (Char *)Malloc(chars*sizeof(Char));
  numszero = (steptr)Malloc(chars*sizeof(steptr));
  numsone = (steptr)Malloc(chars*sizeof(steptr));
  current = (treenumbers)Malloc(spp*sizeof(long));
  order = (treenumbers)Malloc(spp*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  added = (boolean *)Malloc(nonodes*sizeof(boolean));
  ancone = (boolean *)Malloc(chars*sizeof(boolean));
  anczero = (boolean *)Malloc(chars*sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars*sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars*sizeof(boolean));
  wagner = (bitptr)Malloc(words*sizeof(long));
  wagner0 = (bitptr)Malloc(words*sizeof(long));
  zeroanc = (bitptr)Malloc(words*sizeof(long));
  oneanc = (bitptr)Malloc(words*sizeof(long));
  suppsteps = (bitptr)Malloc(words*sizeof(long));
  extras = (steptr)Malloc(chars*sizeof(steptr));
}  /* allocrest */


void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &chars, &nonodes, 1);
  words = chars / bits + 1;
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld characters\n", spp, chars);
  alloctree2(&treenode);
  setuptree2(treenode);
  allocrest();
}  /* doinit */

void inputoptions()
{
  /* input the information on the options */
  Char ch;
  long extranum, i;
  boolean avar, mix;

  if (!firstset)
    samenumsp(&chars, ith);
  extranum = 0;
  avar = false;
  mix = false;
  readoptions(&extranum, "AMW");
  for (i = 0; i < (chars); i++)
    weight[i] = 1;
  for (i = 1; i <= extranum; i++) {
    matchoptions(&ch, "AMW");
    if (ch == 'A') {
      avar = true;
      if (!ancvar) {
	printf("ERROR: ANCESTOR OPTION NOT CHOSEN IN MENU");
	printf(" WITH OPTION %c IN INPUT\n", ch);
	exxit(-1);
      } else
	inputancestors(anczero0, ancone0);
    }
    if (ch == 'M') {
      mix = true;
      if (!mixture) {
	printf("ERROR: MIXTURE OPTION NOT CHOSEN IN MENU");
	printf(" WITH OPTION %c IN INPUT\n", ch);
	exxit(-1);
      } else
	inputmixture(wagner0);
    }
    if (ch == 'W')
      inputweights(chars, weight, &weights);
  }
  if (ancvar && !avar) {
    printf("ERROR: ANCESTOR OPTION CHOSEN IN MENU");
    printf(" WITH NO OPTION A IN INPUT\n");
    exxit(-1);
  }
  if (mixture && !mix) {
    printf("ERROR: MIXTURE OPTION CHOSEN IN MENU");
    printf(" WITH NO OPTION M IN INPUT\n");
    exxit(-1);
  }
  if (weights && printdata)
    printweights(outfile, 0, chars, weight, "Characters");
  for (i = 0; i < (words); i++) {
    if (mixture && mix)
      wagner[i] = wagner0[i];
    else if (allsokal)
      wagner[i] = 0;
    else
      wagner[i] = (1L << (bits + 1)) - (1L << 1);
  }
  if (allsokal)
    fprintf(outfile, "Camin-Sokal parsimony method\n\n");
  if (allwagner)
    fprintf(outfile, "Wagner parsimony method\n\n");
  if (mixture && mix && printdata)
    printmixture(outfile, wagner);
  for (i = 0; i < (chars); i++) {
    if (!ancvar) {
      anczero[i] = true;
      ancone[i] = (((1L << (i % bits + 1)) & wagner[i / bits]) != 0);
    } else {
      anczero[i] = anczero0[i];
      ancone[i] = ancone0[i];
    }
  }
  if (ancvar && avar && printdata)
    printancestors(outfile, anczero, ancone);
  noroot = true;
  questions = false;
  for (i = 0; i < (chars); i++) {
    if (weight[i] > 0) {
      noroot = (noroot && ancone[i] && anczero[i] &&
	  ((((1L << (i % bits + 1)) & wagner[i / bits]) != 0)
            || threshold <= 2.0));
    }
    questions = (questions || (ancone[i] && anczero[i]));
    threshwt[i] = threshold * weight[i];
  }
}  /* inputoptions */


void doinput()
{
  /* reads the input data */
  inputoptions();
  inputdata2(treenode);
}  /* doinput */


void supplement(bitptr suppsteps)
{
  /* determine minimum number of steps more which will
     be added when rest of species are put in tree */
  long i, j, k, l;
  long defone, defzero, a;

  k = 0;
  for (i = 0; i < (words); i++) {
    defone = 0;
    defzero = 0;
    a = 0;
    for (l = 1; l <= bits; l++) {
      k++;
      if (k <= chars) {
	if (!ancone[k - 1])
	  defzero = ((long)defzero) | (1L << l);
	if (!anczero[k - 1])
	  defone = ((long)defone) | (1L << l);
      }
    }
    for (j = 0; j < (spp); j++) {
      defone |= treenode[j]->empstte1[i] & (~treenode[j]->empstte0[i]);
      defzero |= treenode[j]->empstte0[i] & (~treenode[j]->empstte1[i]);
      if (added[j])
	a |= defone & defzero;
    }
    suppsteps[i] = defone & defzero & (~a);
  }
}  /* supplement */


void evaluate(node2 *r)
{
  /* Determines the number of steps needed for a tree.
     This is the minimum number needed to evolve chars on
     this tree */
  long i, stepnum, smaller;
  double sum;

  sum = 0.0;
  for (i = 0; i < (chars); i++) {
    numszero[i] = 0;
    numsone[i] = 0;
  }
  supplement(suppsteps);
  for (i = 0; i < (words); i++)
    zeroanc[i] =fullset;
  full = true;
  postorder(r, fullset, full, wagner, zeroanc);
  cpostorder(r, full, zeroanc, numszero, numsone);
  count(r->fulstte1, zeroanc, numszero, numsone);
  count(suppsteps, zeroanc, numszero, numsone);
  for (i = 0; i < (words); i++)
    zeroanc[i] = 0;
  full = false;
  postorder(r, fullset, full, wagner, zeroanc);
  cpostorder(r, full, zeroanc, numszero, numsone);
  count(r->empstte0, zeroanc, numszero, numsone);
  count(suppsteps, zeroanc, numszero, numsone);
  for (i = 0; i < (chars); i++) {
    smaller = spp * weight[i];
    numsteps[i] = smaller;
    if (anczero[i]) {
      numsteps[i] = numszero[i];
      smaller = numszero[i];
    }
    if (ancone[i] && numsone[i] < smaller)
      numsteps[i] = numsone[i];
    stepnum = numsteps[i] + extras[i];
    if (stepnum <= threshwt[i])
      sum += stepnum;
    else
      sum += threshwt[i];
    guess[i] = '?';
    if (!ancone[i] || (anczero[i] && numszero[i] < numsone[i]))
      guess[i] = '0';
    else if (!anczero[i] || (ancone[i] && numsone[i] < numszero[i]))
      guess[i] = '1';
  }
  if (examined == 0 && mults == 0)
    bestyet = -1.0;
  like = sum;
}  /* evaluate */


void addtraverse(node2 *a, node2 *b, node2 *c, long *m, long *n,
			valptr valyew, placeptr place)
{
  /* traverse all places to add b */
  if (done)
    return;
  if ((*m) <= 2 || !(noroot && (a == root || a == root->next->back))) {
    add3(a, b, c, &root, treenode);
    (*n)++;
    evaluate(root);
    examined++;
    if (examined == howoften) {
      examined = 0;
      mults++;
      if (mults == howmany)
	done = true;
      if (progress) {
	printf("%6ld", mults);
	if (bestyet >= 0)
	  printf("%18.5f", bestyet);
	else
	  printf("         -        ");
	printf("%17ld%20.2f\n", nextree - 1, fracdone * 100);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
    }
    valyew[(*n) - 1] = like;
    place[(*n) - 1] = a->index;
    re_move3(&b, &c, &root, treenode);
  }
  if (!a->tip) {
    addtraverse(a->next->back, b, c, m,n,valyew,place);
    addtraverse(a->next->next->back, b, c, m,n,valyew,place);
  }
}  /* addtraverse */


void addit(long m)
{
  /* adds the species one by one, recursively */
  long n;
  valptr valyew;
  placeptr place;
  long i, j, n1, besttoadd = 0;
  valptr bestval;
  placeptr bestplace;
  double oldfrac, oldfdone, sum, bestsum;

  valyew = (valptr)Malloc(nonodes*sizeof(double));
  bestval = (valptr)Malloc(nonodes*sizeof(double));
  place = (placeptr)Malloc(nonodes*sizeof(long));
  bestplace = (placeptr)Malloc(nonodes*sizeof(long));
  if (simple && !firsttime) {
    n = 0;
    added[order[m - 1] - 1] = true;
    addtraverse(root, treenode[order[m - 1] - 1],
		treenode[spp + m - 2], &m,&n,valyew,place);
    besttoadd = order[m - 1];
    memcpy(bestplace, place, nonodes*sizeof(long));
    memcpy(bestval, valyew, nonodes*sizeof(double));
  } else {
    bestsum = -1.0;
    for (i = 1; i <= (spp); i++) {
      if (!added[i - 1]) {
	n = 0;
	added[i - 1] = true;
	addtraverse(root, treenode[i - 1], treenode[spp + m - 2], &m,
		    &n,valyew,place);
	added[i - 1] = false;
	sum = 0.0;
	for (j = 0; j < (n); j++)
	  sum += valyew[j];
	if (sum > bestsum) {
	  bestsum = sum;
	  besttoadd = i;
	  memcpy(bestplace, place, nonodes*sizeof(long));
	  memcpy(bestval, valyew, nonodes*sizeof(double));
	}
      }
    }
  }
  order[m - 1] = besttoadd;
  memcpy(place, bestplace, nonodes*sizeof(long));
  memcpy(valyew, bestval, nonodes*sizeof(double));
  shellsort(valyew, place, n);
  oldfrac = fracinc;
  oldfdone = fracdone;
  n1 = 0;
  for (i = 0; i < (n); i++) {
    if (valyew[i] <= bestyet || bestyet < 0.0)
      n1++;
  }
  if (n1 > 0)
    fracinc /= n1;
  for (i = 0; i < (n); i++) {
    if (valyew[i] <= bestyet || bestyet < 0.0) {
      current[m - 1] = place[i];
      add3(treenode[place[i] - 1], treenode[besttoadd - 1],
	  treenode[spp + m - 2], &root, treenode);
      added[besttoadd - 1] = true;
      if (m < spp)
	addit(m + 1);
      else {
	if (valyew[i] < bestyet || bestyet < 0.0) {
	  nextree = 1;
	  bestyet = valyew[i];
	}
	if (nextree <= maxtrees) {
	  memcpy(bestorders[nextree - 1], order,
		 spp*sizeof(long));
	  memcpy(bestrees[nextree - 1], current,
		 spp*sizeof(long));
	}
	nextree++;
	firsttime = false;
      }
      re_move3(&treenode[besttoadd - 1], &treenode[spp + m - 2], &root, treenode);
      added[besttoadd - 1] = false;
    }
    fracdone += fracinc;
  }
  fracinc = oldfrac;
  fracdone = oldfdone;
  free(valyew);
  free(bestval);
  free(place);
  free(bestplace);
}  /* addit */


void reroot(node2 *outgroup)
{
  /* reorients tree, putting outgroup in desired position. */
  node2 *p, *q, *newbottom, *oldbottom;

  if (outgroup->back->index == root->index)
    return;
  newbottom = outgroup->back;
  p = treenode[newbottom->index - 1]->back;
  while (p->index != root->index) {
    oldbottom = treenode[p->index - 1];
    treenode[p->index - 1] = p;
    p = oldbottom->back;
  }
  p = root->next;
  q = root->next->next;
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  outgroup->back->back = root->next->next;
  outgroup->back = root->next;
  treenode[newbottom->index - 1] = newbottom;
}  /* reroot */


void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
     each character */

  if (stepbox) {
    putc('\n', outfile);
    writesteps(weights, numsteps);
  }
  if (questions && (!noroot || didreroot))
    guesstates(guess);
  if (ancseq) {
    hypstates(fullset, full, noroot, didreroot, root, wagner,
                zeroanc, oneanc, treenode, guess, garbage);
    putc('\n', outfile);
  }
  if (trout) {
    col = 0;
    treeout2(root, &col, root);
  }
}  /* describe */


void maketree()
{
  /* tree construction recursively by branch and bound */
  long i, j, k;
  node2 *dummy;

  fullset = (1L << (bits + 1)) - (1L << 1);
  if (progress) {
    printf("\nHow many\n");
    printf("trees looked                                       Approximate\n");
    printf("at so far      Length of        How many           percentage\n");
    printf("(multiples     shortest tree    trees this long    searched\n");
    printf("of %4ld):      found so far     found so far       so far\n",
	   howoften);
    printf("----------     ------------     ------------       ------------\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
  }
  done = false;
  mults = 0;
  examined = 0;
  nextree = 1;
  root = treenode[0];
  firsttime = true;
  for (i = 0; i < (spp); i++)
    added[i] = false;
  added[0] = true;
  order[0] = 1;
  k = 2;
  fracdone = 0.0;
  fracinc = 1.0;
  bestyet = -1.0;
  addit(k);
  if (done) {
    if (progress) {
      printf("Search broken off!  Not guaranteed to\n");
      printf(" have found the most parsimonious trees.\n");
    }
    if (treeprint) {
      fprintf(outfile, "Search broken off!  Not guaranteed to\n");
      fprintf(outfile, " have found the most parsimonious\n");
      fprintf(outfile, " trees, but here is what we found:\n");
    }
  }
  if (treeprint) {
    fprintf(outfile, "\nrequires a total of %18.3f\n\n", bestyet);
    if (nextree == 2)
      fprintf(outfile, "One most parsimonious tree found:\n");
    else
      fprintf(outfile, "%5ld trees in all found\n", nextree - 1);
  }
  if (nextree > maxtrees + 1) {
    if (treeprint)
      fprintf(outfile, "here are the first%4ld of them\n", (long)maxtrees);
    nextree = maxtrees + 1;
  }
  if (treeprint)
    putc('\n', outfile);
  for (i = 0; i < (spp); i++)
    added[i] = true;
  for (i = 0; i <= (nextree - 2); i++) {
    for (j = k; j <= (spp); j++)
      add3(treenode[bestrees[i][j - 1] - 1],
	  treenode[bestorders[i][j - 1] - 1], treenode[spp + j - 2],
          &root, treenode);
    if (noroot)
      reroot(treenode[outgrno - 1]);
    didreroot = (outgropt && noroot);
    evaluate(root);
    printree(treeprint, noroot, didreroot, root);
    describe();
    for (j = k - 1; j < (spp); j++)
      re_move3(&treenode[bestorders[i][j] - 1], &dummy, &root, treenode);
  }
  if (progress) {
    printf("\nOutput written to output file\n\n");
    if (trout)
      printf("Trees also written onto tree file\n\n");
  }
  if (ancseq)
    freegarbage(&garbage);
}  /* maketree */


int main(int argc, Char *argv[])
{  /* Penny's branch-and-bound method */
   /* Reads in the number of species, number of characters,
     options and data.  Then finds all most parsimonious trees */
#ifdef MAC
  argc = 1;		/* macsetup("Penny","");		*/
  argv[0] = "Penny";
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
  firstset = true;
  garbage = NULL;
  bits = 8*sizeof(long) - 1;
  doinit();
  if (trout)
    openfile(&outtree,OUTTREE,"output tree file", "w",argv[0],outtreename);
  for (ith = 1; ith <= datasets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      fprintf(outfile, "\nData set # %ld:\n\n",ith);
      if (progress)
        printf("\nData set # %ld:\n",ith);
    }
    maketree();
  }
  FClose(infile);
  FClose(outfile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* Penny's branch-and-bound method */
