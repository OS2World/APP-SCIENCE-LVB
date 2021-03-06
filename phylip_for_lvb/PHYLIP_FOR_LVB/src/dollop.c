
#include "phylip.h"
#include "disc.h"
#include "dollo.h"

/* version 3.6. (c) Copyright 1993-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxtrees        100  /* maximum number of tied trees stored     */

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   inputoptions(void);
void   doinput(void);
void   dollop_count(node *, steptr, steptr);
void   preorder(node *, steptr, steptr, long, boolean, long, bitptr,
	pointptr);
void   evaluate(node *);
void   savetree(void);
void   dollop_addtree(long *);

void   tryadd(node *, node **, node **);
void   addpreorder(node *, node *, node *);
void   tryrearr(node *, node **, boolean *);
void   repreorder(node *, node **, boolean *);
void   rearrange(node **);
void   describe(void);
void   initdollopnode(node **, node **, node *, long, long, long *,
	long *, initops, pointarray, pointarray, Char *, Char *, FILE *);
void   maketree(void);
/* function prototypes */
#endif


extern FILE *infile, *outfile, *intree, *outtree;
extern long spp, chars, nonodes, words, bits, nextree, which;
extern boolean printdata;
extern steptr extras, weight;
extern naym *nayme;

Char infilename[100], outfilename[100], intreename[100], outtreename[100];
node *root;
long col, datasets, ith, j, l, njumble, jumb;
long inseed;
extern boolean ibmpc, ansi;
boolean jumble, usertree, weights, thresh, ancvar, questions, dollo,
               trout,  progress, treeprint, stepbox, ancseq, mulsets,
               firstset;
boolean *ancone, *anczero, *ancone0, *anczero0;
pointptr treenode;   /* pointers to all nodes in tree */
double threshold;
double *threshwt;
longer seed;
long *enterorder;
double **fsteps;
steptr numsteps;
bestelm *bestrees;
Char *guess;
gbit *garbage;
Char progname[50];

/* Variables for treeread */
boolean goteof, firsttree, haslengths, phirst;
pointarray nodep;
node *grbg;
long *zeros;

/* Local variables for maketree, propagated globally for C version: */
long minwhich;
double like, bestyet, bestlike, bstlike2, minsteps;
boolean lastrearr;
double nsteps[maxuser];
node *there;
long fullset;
bitptr zeroanc, oneanc;
long *place;
Char ch;
boolean *names;
steptr numsone, numszero;
bitptr steps;


void getoptions()
{
  /* interactively set options */
  long inseed0;
  Char ch;

  fprintf(outfile,"\nDollo and polymorphism parsimony algorithm,");
  fprintf(outfile," version %s\n\n",VERSION);
  putchar('\n');
  ancvar = false;
  dollo = true;
  jumble = false;
  njumble = 1;
  thresh = false;
  threshold = spp;
  trout = true;
  usertree = false;
  goteof = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  for (;;) {
    printf(ansi ? "\033[2J\033[H" : "\n");
    printf("\nDollo and polymorphism parsimony algorithm, version %s\n\n",
           VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
	   (usertree ? "No, use user trees in input file" : "Yes"));
    printf("  P                     Parsimony method?  %s\n",
	   dollo ? "Dollo" : "Polymorphism");
    if (!usertree) {
      printf("  J     Randomize input order of species?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
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
    printf("\nAre these settings correct? ");
    printf("(type Y or the letter for one to change)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("APJTUM1234560",ch) != NULL){
      switch (ch) {
	
      case 'A':
	ancvar = !ancvar;
	break;
	
      case 'P':
	dollo = !dollo;
	break;
	
      case 'J':
	jumble = !jumble;
	if (jumble)
	  initjumble(&inseed, &inseed0, seed, &njumble);
	else njumble = 1;
	break;
	
      case 'T':
	thresh = !thresh;
	if (thresh)
	  initthreshold(&threshold);
	break;
	
      case 'U':
	usertree = !usertree;
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
}  /* getoptions */


void allocrest()
{
  long i;

  extras = (steptr)Malloc(chars*sizeof(long));
  weight = (steptr)Malloc(chars*sizeof(long));
  threshwt = (double *)Malloc(chars*sizeof(double));
  if (usertree) {
    fsteps = (double **)Malloc(maxuser*sizeof(double *));
    for (i = 1; i <= maxuser; i++)
      fsteps[i - 1] = (double *)Malloc(chars*sizeof(double));
  }
  bestrees = (bestelm *) Malloc(maxtrees*sizeof(bestelm));
  for (i = 1; i <= maxtrees; i++)
    bestrees[i - 1].btree = (long *)Malloc(nonodes*sizeof(long));
  numsteps = (steptr)Malloc(chars*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  enterorder = (long *)Malloc(spp*sizeof(long));
  place = (long *)Malloc(nonodes*sizeof(long));
  ancone = (boolean *)Malloc(chars*sizeof(boolean));
  anczero = (boolean *)Malloc(chars*sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars*sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars*sizeof(boolean));
  numsone = (steptr)Malloc(chars*sizeof(long));
  numszero = (steptr)Malloc(chars*sizeof(long));
  guess = (Char *)Malloc(chars*sizeof(Char));
  zeroanc = (bitptr)Malloc(words*sizeof(long));
  oneanc = (bitptr)Malloc(words*sizeof(long));
}  /* allocrest */


void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &chars, &nonodes, 1);
  words = chars / bits + 1;
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  characters\n\n", spp, chars);
  alloctree(&treenode);
  setuptree(treenode);
  allocrest();
}  /* doinit */


void inputoptions()
{
  /* input the information on the options */
  Char ch;
  long extranum, i;
  boolean avar;

  if (!firstset)
    samenumsp(&chars, ith);
  extranum = 0;
  avar = false;
  readoptions(&extranum, "AW");
  for (i = 0; i < (chars); i++)
    weight[i] = 1;
  for (i = 1; i <= extranum; i++) {
    matchoptions(&ch, "AW");
    if (ch == 'A') {
      avar = true;
      if (!ancvar) {
	printf("ERROR: ANCESTOR OPTION NOT CHOSEN IN MENU");
	printf(" WITH OPTION %c IN INPUT\n", ch);
	exxit(-1);
      } else
	inputancestors(anczero0, ancone0);
    }
    if (ch == 'W')
      inputweights(chars, weight, &weights);
  }
  if (ancvar && !avar) {
    puts("ERROR: ANCESTOR OPTION CHOSEN IN MENU WITH NO OPTION A IN INPUT");
    exxit(-1);
  }
  if (dollo)
    fprintf(outfile, "Dollo");
  else
    fprintf(outfile, "Polymorphism");
  fprintf(outfile, " parsimony method\n\n");
  if (weights && printdata)
    printweights(outfile, 0, chars, weight, "Characters");
  for (i = 0; i < (chars); i++) {
    if (!ancvar) {
      anczero[i] = true;
      ancone[i] = false;
    } else {
      anczero[i] = anczero0[i];
      ancone[i] = ancone0[i];
    }
  }
  if (ancvar && avar && printdata)
    printancestors(outfile, anczero, ancone);
  questions = false;
  for (i = 0; i < (chars); i++) {
    questions = (questions || (ancone[i] && anczero[i]));
    threshwt[i] = threshold * weight[i];
  }
}  /* inputoptions */


void doinput()
{
  /* reads the input data */
  inputoptions();
  inputdata(treenode, dollo, printdata, outfile);
}  /* doinput */


void dollop_count(node *p, steptr numsone, steptr numszero)
{
  /* counts the number of steps in a fork of the tree.
     The program spends much of its time in this PROCEDURE */
  long i, j, l;

  if (dollo) {
    for (i = 0; i < (words); i++)
      steps[i] = (treenode[p->back->index - 1]->stateone[i] &
                  p->statezero[i] & zeroanc[i]) |
          (treenode[p->back->index - 1]->statezero[i] & p->stateone[i] &
           fullset & (~zeroanc[i]));
  } else {
    for (i = 0; i < (words); i++)
      steps[i] = treenode[p->back->index - 1]->stateone[i] &
                 treenode[p->back->index - 1]->statezero[i] & p->stateone[i] &
                 p->statezero[i];
  }
  j = 1;
  l = 0;
  for (i = 0; i < (chars); i++) {
    l++;
    if (l > bits) {
      l = 1;
      j++;
    }
    if (((1L << l) & steps[j - 1]) != 0) {
      if (((1L << l) & zeroanc[j - 1]) != 0)
        numszero[i] += weight[i];
      else
        numsone[i] += weight[i];
    }
  }
}  /* dollop_count */


void preorder(node *p, steptr numsone, steptr numszero, long words,
		boolean dollo, long fullset, bitptr zeroanc, pointptr treenode)
{
  /* go back up tree setting up and counting interior node
     states */

  if (!p->tip) {
    correct(p, fullset, dollo, zeroanc, treenode);
    preorder(p->next->back, numsone,numszero, words, dollo, fullset,
               zeroanc, treenode);
    preorder(p->next->next->back, numsone,numszero, words, dollo, fullset,
               zeroanc, treenode);
  }
  if (p->back != NULL)
    dollop_count(p, numsone,numszero);
}  /* preorder */


void evaluate(node *r)
{
  /* Determines the number of losses or polymorphisms needed
     for a tree. This is the minimum number needed to evolve
     chars on this tree */
  long i, stepnum, smaller;
  double sum, term;

  sum = 0.0;
  for (i = 0; i < (chars); i++) {
    numszero[i] = 0;
    numsone[i] = 0;
  }
  for (i = 0; i < (words); i++)
    zeroanc[i] = fullset;
  postorder(r);
  preorder(r, numsone, numszero, words, dollo, fullset, zeroanc, treenode);
  for (i = 0; i < (words); i++)
    zeroanc[i] = 0;
  postorder(r);
  preorder(r, numsone, numszero, words, dollo, fullset, zeroanc, treenode);
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
      term = stepnum;
    else
      term = threshwt[i];
    sum += term;
    if (usertree && which <= maxuser)
      fsteps[which - 1][i] = term;
    guess[i] = '?';
    if (!ancone[i] || (anczero[i] && numszero[i] < numsone[i]))
      guess[i] = '0';
    else if (!anczero[i] || (ancone[i] && numsone[i] < numszero[i]))
      guess[i] = '1';
  }
  if (usertree && which <= maxuser) {
    nsteps[which - 1] = sum;
    if (which == 1) {
      minwhich = 1;
      minsteps = sum;
    } else if (sum < minsteps) {
      minwhich = which;
      minsteps = sum;
    }
  }
  like = -sum;
}  /* evaluate */


void savetree()
{
  /* record in place where each species has to be
     added to reconstruct this tree */
  long i, j;
  node *p;
  boolean done;

  for (i = 0; i < (nonodes); i++)
    place[i] = 0;
  place[root->index - 1] = 1;
  for (i = 1; i <= (spp); i++) {
    p = treenode[i - 1];
    while (place[p->index - 1] == 0) {
      place[p->index - 1] = i;
      p = p->back;
      if (p != NULL)
        p = treenode[p->index - 1];
    }
    if (i > 1) {
      place[i - 1] = place[p->index - 1];
      j = place[p->index - 1];
      done = false;
      while (!done) {
        place[p->index - 1] = spp + i - 1;
        p = treenode[p->index - 1];
        p = p->back;
        done = (p == NULL);
        if (!done)
          done = (place[p->index - 1] != j);
      }
    }
  }
}  /* savetree */


void dollop_addtree(long *pos)
{
  /*puts tree from ARRAY place in its proper position
    in ARRAY bestrees */
  long i;
  for (i =nextree - 1; i >= (*pos); i--) {
    memcpy(bestrees[i].btree, bestrees[i - 1].btree, spp*sizeof(long));
    bestrees[i].gloreange = bestrees[i - 1].gloreange;
    bestrees[i].locreange = bestrees[i - 1].locreange;
    bestrees[i].collapse = bestrees[i - 1].collapse;
  }
  for (i = 0; i < (spp); i++)
    bestrees[(*pos) - 1].btree[i] = place[i];
  nextree++;
}  /* dollop_addtree */


void tryadd(node *p, node **item, node **nufork)
{
  /* temporarily adds one fork and one tip to the tree.
     if the location where they are added yields greater
     "likelihood" than other locations tested up to that
     time, then keeps that location as there */
  long pos;
  boolean found;

  add(p, *item, *nufork, &root, treenode);
  evaluate(root);
  if (lastrearr) {
    if (like >= bstlike2) {
      savetree();
      if (like > bstlike2) {
        bestlike = bstlike2 = like;
        pos = 1;
        nextree = 1;
        dollop_addtree(&pos);
      } else {
        pos = 0;
        findtree(&found, &pos, nextree, place, bestrees);
		/* findtree calls for a bestelm* but is getting	*/
		/* a long**, LM					*/
        if (!found) {
          if (nextree <= maxtrees)
            dollop_addtree(&pos);
        }
      }
    }
  }
  if (like > bestyet) {
    bestyet = like;
    there = p;
  }
  re_move(item, nufork, &root, treenode);
}  /* tryadd */


void addpreorder(node *p, node *item_, node *nufork_)
{
  /* traverses a binary tree, calling PROCEDURE tryadd
     at a node before calling tryadd at its descendants */
  node *item= item_;
  node *nufork = nufork_;

  if (p == NULL)
    return;
  tryadd(p, &item,&nufork);
  if (!p->tip) {
    addpreorder(p->next->back, item, nufork);
    addpreorder(p->next->next->back, item, nufork);
  }
}  /* addpreorder */


void tryrearr(node *p, node **r, boolean *success)
{
  /* evaluates one rearrangement of the tree.
     if the new tree has greater "likelihood" than the old
     one sets success := TRUE and keeps the new tree.
     otherwise, restores the old tree */
  node *frombelow, *whereto, *forknode;
  double oldlike;

  if (p->back == NULL)
    return;
  forknode = treenode[p->back->index - 1];
  if (forknode->back == NULL)
    return;
  oldlike = bestyet;
  if (p->back->next->next == forknode)
    frombelow = forknode->next->next->back;
  else
    frombelow = forknode->next->back;
  whereto = forknode->back;
  re_move(&p, &forknode, &root, treenode);
  add(whereto, p, forknode, &root, treenode);
  evaluate(*r);
  if (like <= oldlike) {
    re_move(&p, &forknode, &root, treenode);
    add(frombelow, p, forknode, &root, treenode);
  } else {
    (*success) = true;
    bestyet = like;
  }
}  /* tryrearr */


void repreorder(node *p, node **r, boolean *success)
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
     at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p, r,success);
  if (!p->tip) {
    repreorder(p->next->back, r,success);
    repreorder(p->next->next->back, r,success);
  }
}  /* repreorder */


void rearrange(node **r_)
{
  /* traverses the tree (preorder), finding any local
     rearrangement which decreases the number of steps.
     if traversal succeeds in increasing the tree's
     "likelihood", PROCEDURE rearrange runs traversal again */
  node **r         = r_;
  boolean success  = true;

  while (success) {
    success = false;
    repreorder(*r, r,&success);
  }
}  /* rearrange */


void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
     each character */

  if (treeprint)
    fprintf(outfile, "\nrequires a total of %10.3f\n", -like);
  if (stepbox) {
    putc('\n', outfile);
    writesteps(weights, dollo, numsteps);
  }
  if (questions)
    guesstates(guess);
  if (ancseq) {
    hypstates(fullset, dollo, guess, treenode, root, garbage, zeroanc, oneanc);
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (trout) {
    col = 0;
    treeout(root, nextree, &col, root);
  }
}  /* describe */


void initdollopnode(node **p, node **grbg, node *q, long len, long nodei,
			long *ntips, long *parens, initops whichinit,
			pointarray treenode, pointarray nodep, Char *str, Char *ch,
			FILE *intree)
{
  /* initializes a node */
  /* LM 7/27  I added this function and the commented lines around */
  /* treeread() to get the program running, but all 4 move programs*/
  /* are improperly integrated into the v4.0 support files.  As is */
  /* this is a patchwork function         */
  switch (whichinit) {
  case bottom:
    gnutreenode(grbg, p, nodei, chars, zeros);
    treenode[nodei - 1] = *p;
    break;
  case nonbottom:
    gnutreenode(grbg, p, nodei, chars, zeros);
    break;
  case tip:
    match_names_to_data (str, treenode, p, spp);
    break;
  default:      /*cases hslength,hsnolength,treewt,unittrwt,iter,*/
    break;      /*length should never occur                      */
  }
} /* initdollopnode */


void maketree()
{
  /* constructs a binary tree from the pointers in treenode.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, numtrees, nextnode;
  double gotlike;
  node *item, *nufork, *dummy, *p;

  steps = (bitptr)Malloc(words*sizeof(long));
  fullset = (1L << (bits + 1)) - (1L << 1);
  if (!usertree) {
    for (i = 1; i <= (spp); i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    root = treenode[enterorder[0] - 1];
    add(treenode[enterorder[0] - 1], treenode[enterorder[1] - 1],
        treenode[spp], &root, treenode);
    if (progress) {
      printf("Adding species:\n");
      writename(0, 2, enterorder);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    lastrearr = false;
    for (i = 3; i <= (spp); i++) {
      bestyet = -350.0 * spp * chars;
      item = treenode[enterorder[i - 1] - 1];
      nufork = treenode[spp + i - 2];
      addpreorder(root, item, nufork);
      add(there, item, nufork, &root, treenode);
      like = bestyet;
      rearrange(&root);
      if (progress) {
        writename(i - 1, 1, enterorder);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
      lastrearr = (i == spp);
      if (lastrearr) {
        if (progress) {
          printf("\nDoing global rearrangements\n");
          printf("  !");
          for (j = 1; j <= (nonodes); j++)
            putchar('-');
          printf("!\n");
#ifdef WIN32
          phyFillScreenColor();
#endif
        }
        bestlike = bestyet;
        if (jumb == 1) {
          bstlike2 = bestlike;
          nextree = 1;
        }
        do {
          if (progress)
            printf("   ");
          gotlike = bestlike;
          for (j = 0; j < (nonodes); j++) {
            bestyet = - 350.0 * spp * chars;
            item = treenode[j];
            if (item != root) {
              nufork = treenode[j]->back;
              re_move(&item, &nufork, &root, treenode);
              there = root;
              addpreorder(root, item, nufork);
              add(there, item, nufork, &root, treenode);
            }
            if (progress) {
              putchar('.');
              fflush(stdout);
            }
          }
          if (progress) {
            putchar('\n');
#ifdef WIN32
            phyFillScreenColor();
#endif

          }
        } while (bestlike > gotlike);
      }
    }
    if (progress)
      putchar('\n');
    for (i = spp - 1; i >= 1; i--)
      re_move(&treenode[i], &dummy, &root, treenode);
    if (jumb == njumble) {
      if (treeprint) {
        putc('\n', outfile);
        if (nextree == 2)
          fprintf(outfile, "One most parsimonious tree found:\n");
        else
          fprintf(outfile, "%6ld trees in all found\n", nextree - 1);
      }
      if (nextree > maxtrees + 1) {
        if (treeprint)
          fprintf(outfile, "here are the first%4ld of them\n", (long)maxtrees);
        nextree = maxtrees + 1;
      }
      if (treeprint)
        putc('\n', outfile);
      for (i = 0; i <= (nextree - 2); i++) {
        root = treenode[0];
        add(treenode[0], treenode[1], treenode[spp], &root, treenode);
        for (j = 3; j <= spp; j++) {
          add(treenode[bestrees[i].btree[j - 1] - 1], treenode[j - 1],
              treenode[spp + j - 2], &root, treenode);}
        evaluate(root);
        printree(1.0, treeprint, root);
        describe();
        for (j = 1; j < (spp); j++)
          re_move(&treenode[j], &dummy, &root, treenode);
      }
    }
  } else {
    openfile(&intree,INTREE,"input tree file", "r",progname,intreename);
    fscanf(intree, "%ld%*[^\n]", &numtrees);
    getc(intree);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n");
    }
    names = (boolean *)Malloc(spp*sizeof(boolean));
    which = 1;
    firsttree = true;                       /**/
    nodep = NULL;                           /**/
    nextnode = 0;                           /**/
    haslengths = 0;                         /**/
    phirst = 0;				    /**/
    zeros = (long *)Malloc(chars*sizeof(long));         /**/
    for (i = 0; i < chars; i++)             /**/
      zeros[i] = 0;                         /**/
    while (which <= numtrees) {
      treeread(intree, &root, treenode, &goteof, &firsttree,
                nodep, intreename, "Dollop", &nextnode, &haslengths,
                &grbg, initdollopnode); /*debug*/
      for (i = spp; i < (nonodes); i++) {
        p = treenode[i];
        for (j = 1; j <= 3; j++) {
          p->stateone = (bitptr)Malloc(words*sizeof(long));
          p->statezero = (bitptr)Malloc(words*sizeof(long));
          p = p->next;
        }
      } /* debug: see comment at initdollopnode() */
      if (treeprint)
        fprintf(outfile, "\n\n");
      evaluate(root);
      printree(1.0, treeprint, root);
      describe();
      which++;
    }
    FClose(intree);
    fprintf(outfile, "\n\n");
    if (numtrees > 1 && chars > 1)
      standev(numtrees, which, minwhich, minsteps, nsteps, fsteps);
    free(names);
  }
  if (jumb == njumble) {
    if (progress) {
      printf("Output written to output file\n\n");
      if (trout)
        printf("Trees also written onto tree file\n\n");
    }
    free(steps);
    if (ancseq)
      freegarbage(&garbage);
  }
}  /* maketree */


int main(int argc, Char *argv[])
{  /* Dollo or polymorphism parsimony by uphill search */
#ifdef MAC
  argc = 1;		/* macsetup("Dollop","");		*/
  argv[0] = "Dollop";
#endif
#ifdef WIN32
  phySetConsoleAttributes();
  phyClearScreen();
#endif
  /* reads in spp, chars, and the data. Then calls maketree to
     construct the tree */
  strcpy(progname,argv[0]);
  openfile(&infile,INFILE,"input file", "r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file", "w",argv[0],outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  garbage = NULL;
  mulsets = false;
  datasets = 1;
  firstset = true;
  bits = 8*sizeof(long) - 1;
  doinit();
  if (trout)
    openfile(&outtree,OUTTREE,"output tree file", "w",argv[0],outtreename);
  for (ith = 1; ith <= (datasets); ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n\n",ith);
      if (progress)
        printf("\nData set # %ld:\n",ith);
    }
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
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* Dollo or polymorphism parsimony by uphill search */
