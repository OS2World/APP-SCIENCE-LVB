#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1993-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxtrees        100   /* maximum number of tied trees stored     */

typedef boolean *boolptr;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   initdnacompnode(node **, node **, node *, long, long, long *,
		long *, initops, pointarray, pointarray, Char *, Char *, FILE *);
void   inputoptions(void);
void   makeweights(void);
void   doinput(void);
void   mincomp(long );
void   evaluate(node *);
void   localsavetree(void);

void   tryadd(node *, node *, node *);
void   addpreorder(node *, node *, node *);
void   tryrearr(node *, boolean *);
void   repreorder(node *, boolean *);
void   rearrange(node **);
void   describe(void);
void   initboolnames(node *, boolean *);
void   maketree(void);
void   freerest(void);
/* function prototypes */
#endif


extern FILE *infile, *outfile, *intree, *outtree;
extern long spp, nonodes, endsite, outgrno, nextree, which;
extern boolean ibmpc, ansi, interleaved, printdata, outgropt, treeprint;
extern sequence y;
extern steptr weight, alias, location, ally;
extern naym *nayme;

Char infilename[100], outfilename[100], intreename[100], outtreename[100];
node *root, *p;
long chars, col, datasets, ith, njumble, jumb;
long inseed;
boolean jumble, usertree, trout, weights,
               progress, stepbox, ancseq, firstset, mulsets;
steptr oldweight, necsteps;
pointarray treenode;   /* pointers to all nodes in tree */
long *enterorder;
Char basechar[32]="ACMGRSVTWYHKDBNO???????????????";
bestelm *bestrees;
boolean dummy;
longer seed;
gbases *garbage;
Char ch;
Char progname[20];
long *zeros;

/* Local variables for maketree, propogated globally for C version: */
  long maxwhich;
  double like, maxsteps, bestyet, bestlike, bstlike2;
  boolean lastrearr, recompute;
  double nsteps[maxuser];
  long **fsteps;
  node *there;
  long *place;
  boolptr in_tree;
  baseptr nothing;
  node *temp, *temp1;
  node *grbg;

void getoptions()
{
  /* interactively set options */
  long inseed0;
  Char ch;

  fprintf(outfile, "\nDNA compatibility algorithm, version %s\n\n",VERSION);
  putchar('\n');
  jumble = false;
  njumble = 1;
  outgrno = 1;
  outgropt = false;
  trout = true;
  usertree = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  interleaved = true;
  for (;;) {
    if (ansi)
      printf("\033[2J\033[H");
    else
      putchar('\n');
    printf("\nDNA compatibility algorithm, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
	   (usertree ? "No, use user trees in input file" : "Yes"));
    if (!usertree) {
      printf("  J   Randomize input order of sequences?");
      if (jumble) {
        printf(
         "  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      }
      else
        printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at sequence number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
	   (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
	   ibmpc ? "IBM PC" : ansi  ? "ANSI"   : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
	   (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
	   (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
	   (treeprint ? "Yes" : "No"));
    printf("  4  Print steps & compatibility at sites  %s\n",
	   (stepbox ? "Yes" : "No"));
    printf("  5  Print sequences at all nodes of tree  %s\n",
	   (ancseq ? "Yes" : "No"));
    printf("  6       Write out trees onto tree file?  %s\n",
	   (trout ? "Yes" : "No"));
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("JOTUMI1234560",ch) != NULL){
      switch (ch) {
	
      case 'J':
	jumble = !jumble;
	if (jumble)
	  initjumble(&inseed, &inseed0, seed, &njumble);
	else njumble = 1;
	break;
	
      case 'O':
	outgropt = !outgropt;
	if (outgropt)
	  initoutgroup(&outgrno, spp);
	break;
	
      case 'U':
	usertree = !usertree;
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

  y = (Char **)Malloc(spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    y[i] = (Char *)Malloc(chars*sizeof(Char));
  bestrees = (bestelm *) Malloc(maxtrees*sizeof(bestelm));
  for (i = 1; i <= maxtrees; i++)
    bestrees[i - 1].btree = (long *)Malloc(nonodes*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  weight = (steptr)Malloc(chars*sizeof(long));
  oldweight = (steptr)Malloc(chars*sizeof(long));
  enterorder = (long *)Malloc(spp*sizeof(long));
  necsteps = (steptr)Malloc(chars*sizeof(long));
  alias = (steptr)Malloc(chars*sizeof(long));
  ally = (steptr)Malloc(chars*sizeof(long));
  location = (steptr)Malloc(chars*sizeof(long));
  place = (long *)Malloc((2*spp-1)*sizeof(long));
  in_tree = (boolptr)Malloc(chars*sizeof(boolean));
}  /* allocrest */


void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &chars, &nonodes, 1);
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, chars);
  alloctree(&treenode, nonodes, usertree);
  allocrest();
}  /* doinit */


void initdnacompnode(node **p, node **grbg, node *q, long len, long nodei,
			long *ntips, long *parens, initops whichinit,
			pointarray treenode, pointarray nodep, Char *str, Char *ch,
			FILE *intree)
{
  /* initializes a node */

  switch (whichinit) {
  case bottom:
    gnutreenode(grbg, p, nodei, endsite, zeros);
    treenode[nodei - 1] = *p;
    break;
  case nonbottom:
    gnutreenode(grbg, p, nodei, endsite, zeros);
    break;
  case tip:
    match_names_to_data (str, treenode, p, spp);
    break;
  default:
    break;
  }
} /* initdnacompnode */


void inputoptions()
{
  /* input the information on the options */
  long i;

  if (!firstset)
    samenumsp(&chars, ith);
  for (i = 0; i < chars; i++)
    weight[i] = 1;
  if (weights) {
    inputweights(chars, weight, &weights);
    printweights(outfile, 0, chars, weight, "Sites");
  }
}  /* inputoptions */


void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= chars; i++) {
    alias[i - 1] = i;
    oldweight[i - 1] = weight[i - 1];
    ally[i - 1] = i;
  }
  sitesort(chars, weight);
  sitecombine(chars);
  sitescrunch(chars);
  endsite = 0;
  for (i = 1; i <= chars; i++) {
    if (ally[i - 1] == i)
      endsite++;
  }
  for (i = 1; i <= endsite; i++)
    location[alias[i - 1] - 1] = i;
  zeros = (long *)Malloc(endsite*sizeof(long));
  for (i = 0; i < endsite; i++)
    zeros[i] = 0;
}  /* makeweights */


void doinput()
{
  /* reads the input data */
  inputoptions();
  inputdata(chars);
  makeweights();
  makevalues(treenode, zeros, usertree);
  allocnode(&temp, zeros, endsite);
  allocnode(&temp1, zeros, endsite);
}  /* doinput */


void mincomp(long n)
{
  /* computes for each site the minimum number of steps
    necessary to accomodate those species already
    in the analysis, adding in species n */
  long i, j, k, l, m;
  bases b;
  long s;
  boolean allowable, deleted;

  in_tree[n - 1] = true;
  for (i = 0; i < endsite; i++)
    necsteps[i] = 3;
  for (m = 0; m <= 31; m++) {
    s = 0;
    l = -1;
    k = m;
    for (b = A; (long)b <= (long)O; b = (bases)((long)b + 1)) {
      if ((k & 1) == 1) {
        s |= 1L << ((long)b);
        l++;
      }
      k /= 2;
    }
    for (j = 0; j < endsite; j++) {
      allowable = true;
      i = 1;
      while (allowable && i <= spp) {
        if (in_tree[i - 1] && treenode[i - 1]->base[j] != 0) {
          if ((treenode[i - 1]->base[j] & s) == 0)
            allowable = false;
        }
        i++;
      }
      if (allowable) {
        if (l < necsteps[j])
          necsteps[j] = l;
      }
    }
  }
  for (j = 0; j < endsite; j++) {
    deleted = false;
    for (i = 0; i < spp; i++) {
      if (in_tree[i] && treenode[i]->base[j] == 0)
        deleted = true;
    }
    if (deleted)
      necsteps[j]++;
  }
  for (i = 0; i < endsite; i++)
    necsteps[i] *= weight[i];
}  /* mincomp */


void evaluate(node *r)
{
  /* determines the number of steps needed for a tree. this is
    the minimum number of steps needed to evolve sequences on
    this tree */
  long i, term;
  double sum;

   sum = 0.0;
   for (i = 0; i < endsite; i++) {
    if (r->numsteps[i] == necsteps[i])
      term = weight[i];
    else
      term = 0;
    sum += term;
    if (usertree && which <= maxuser)
      fsteps[which - 1][i] = term;
  }
  if (usertree && which <= maxuser) {
    nsteps[which - 1] = sum;
    if (which == 1) {
      maxwhich = 1;
      maxsteps = sum;
    } else if (sum > maxsteps) {
      maxwhich = which;
      maxsteps = sum;
    }
  }
  like = sum;
}  /* evaluate */


void localsavetree()
{
  /* record in place where each species has to be
    added to reconstruct this tree */
  long i, j;
  node *p;
  boolean done;

  reroot(treenode[outgrno - 1], root);
  savetraverse(root);
  for (i = 0; i < nonodes; i++)
    place[i] = 0;
  place[root->index - 1] = 1;
  for (i = 1; i <= spp; i++) {
    p = treenode[i - 1];
    while (place[p->index - 1] == 0) {
      place[p->index - 1] = i;
      while (!p->bottom)
        p = p->next;
      p = p->back;
    }
    if (i > 1) {
      place[i - 1] = place[p->index - 1];
      j = place[p->index - 1];
      done = false;
      while (!done) {
        place[p->index - 1] = spp + i - 1;
        while (!p->bottom)
          p = p->next;
        p = p->back;
        done = (p == NULL);
        if (!done)
          done = (place[p->index - 1] != j);
      }
    }
  }
}  /* localsavetree */


void tryadd(node *p, node *item, node *nufork)
{
  /* temporarily adds one fork and one tip to the tree.
    if the location where they are added yields greater
    "likelihood" than other locations tested up to that
    time, then keeps that location as there */
  long pos;
  boolean found;
  node *rute, *q;

  if (p == root)
    fillin(temp, item, p);
  else {
    fillin(temp1, item, p);
    fillin(temp, temp1, p->back);
  }
  evaluate(temp);
  if (lastrearr) {
    if (like < bestlike) {
      if (item == nufork->next->next->back) {
        q = nufork->next;
        nufork->next = nufork->next->next;
        nufork->next->next = q;
        q->next = nufork;
      }
    } else if (like >= bstlike2) {
      recompute = false;
      add(p, item, nufork, &root, recompute, treenode, &grbg, zeros);
      rute = root->next->back;
      localsavetree();
      reroot(rute, root);
      if (like > bstlike2) {
        bestlike = bstlike2 = like;
        pos = 1;
        nextree = 1;
        addtree(pos, &nextree, dummy, place, bestrees);
      } else {
        pos = 0;
        findtree(&found, &pos, nextree, place, bestrees);
        if (!found) {
          if (nextree <= maxtrees)
            addtree(pos, &nextree, dummy, place, bestrees);
        }
      }
      re_move(item, &nufork, &root, recompute, treenode, &grbg, zeros);
      recompute = true;
    }
  }
  if (like > bestyet) {
    bestyet = like;
    there = p;
  }
}  /* tryadd */


void addpreorder(node *p, node *item, node *nufork)
{
  /* traverses a binary tree, calling PROCEDURE tryadd
    at a node before calling tryadd at its descendants */

  if (p == NULL)
    return;
  tryadd(p, item, nufork);
  if (!p->tip) {
    addpreorder(p->next->back, item, nufork);
    addpreorder(p->next->next->back, item, nufork);
  }
}  /* addpreorder */


void tryrearr(node *p, boolean *success)
{
  /* evaluates one rearrangement of the tree.
    if the new tree has greater "likelihood" than the old
    one sets success := TRUE and keeps the new tree.
    otherwise, restores the old tree */
  node *frombelow, *whereto, *forknode, *q;
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
  whereto = treenode[forknode->back->index - 1];
  if (whereto->next->back == forknode)
    q = whereto->next->next->back;
  else
    q = whereto->next->back;
  fillin(temp1, frombelow, q);
  fillin(temp, temp1, p);
  fillin(temp1, temp, whereto->back);
  evaluate(temp1);
  if (like <= oldlike) {
    if (p != forknode->next->next->back)
      return;
    q = forknode->next;
    forknode->next = forknode->next->next;
    forknode->next->next = q;
    q->next = forknode;
    return;
  }
  recompute = false;
  re_move(p, &forknode, &root, recompute, treenode, &grbg, zeros);
  fillin(whereto, whereto->next->back, whereto->next->next->back);
  recompute = true;
  add(whereto, p, forknode, &root, recompute, treenode, &grbg, zeros);
  *success = true;
  bestyet = like;
}  /* tryrearr */


void repreorder(node *p, boolean *success)
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
    at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p,success);
  if (!p->tip) {
    repreorder(p->next->back,success);
    repreorder(p->next->next->back,success);
  }
}  /* repreorder */


void rearrange(node **r)
{
  /* traverses the tree (preorder), finding any local
    rearrangement which decreases the number of steps.
    if traversal succeeds in increasing the tree's
    "likelihood", PROCEDURE rearrange runs traversal again */
  boolean success=true;

  while (success) {
    success = false;
    repreorder(*r,&success);
  }
}  /* rearrange */


void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
    each site and table of compatibilities */
  long i, j, k;

  if (treeprint) {
    fprintf(outfile, "\ntotal number of compatible sites is ");
    fprintf(outfile, "%10.1f\n", like);
  }
  if (stepbox) {
    writesteps(chars, weights, oldweight, root);
    fprintf(outfile,
            "\n compatibility (Y or N) of each site with this tree:\n\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%ld", i);
    fprintf(outfile, "\n     *----------\n");
    for (i = 0; i <= (chars / 10); i++) {
      putc(' ', outfile);
      fprintf(outfile, "%3ld !", i * 10);
      for (j = 0; j <= 9; j++) {
        k = i * 10 + j;
        if (k > 0 && k <= chars) {
          if (root->numsteps[location[ally[k - 1] - 1] - 1] ==
              necsteps[location[ally[k - 1] - 1] - 1]) {
            if (oldweight[k - 1] > 0)
              putc('Y', outfile);
            else
              putc('y', outfile);
          } else {
            if (oldweight[k - 1] > 0)
              putc('N', outfile);
            else
              putc('n', outfile);
          }
        } else
          putc(' ', outfile);
      }
      putc('\n', outfile);
    }
  }
  if (ancseq) {
    hypstates(chars, root, treenode, &garbage, basechar);
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (trout) {
    col = 0;
    treeout(root, nextree, &col, root);
  }
}  /* describe */


void initboolnames(node *p, boolean *names)
{
  /* sets BOOLEANs that indicate tips */
  node *q;

  if (p->tip) {
    names[p->index - 1] = true;
    return;
  }
  q = p->next;
  while (q != p) {
    initboolnames(q->back, names);
    q = q->next;
  }
}  /* initboolnames */


void maketree()
{
  /* constructs a binary tree from the pointers in treenode.
    adds each node at location which yields highest "likelihood"
    then rearranges the tree for greatest "likelihood" */
  long i, j, numtrees, num, nextnode;
  boolean firsttree, goteof, haslengths;
  double gotlike, wt, sumw, sum, sum2, sd;
  node *item, *nufork, *dummy;
  long TEMP2;
  pointarray nodep;
  boolean *names;

  if (!usertree) {
    recompute = true;
    for (i = 0; i < spp; i++)
      in_tree[i] = false;
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    root = treenode[enterorder[0] - 1];
    add(treenode[enterorder[0] - 1], treenode[enterorder[1] - 1],
        treenode[spp], &root, recompute, treenode, &grbg, zeros);
    if (progress) {
      printf("Adding species:\n");
      writename(0, 2, enterorder);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    in_tree[0] = true;
    in_tree[1] = true;
    lastrearr = false;
    for (i = 3; i <= spp; i++) {
      mincomp(i);
      bestyet = -350.0 * spp * chars;
      item = treenode[enterorder[i - 1] - 1];
      nufork = treenode[spp + i - 2];
      there = root;
      addpreorder(root, item, nufork);
      add(there, item, nufork, &root, recompute, treenode, &grbg, zeros);
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
          for (j = 1; j <= nonodes; j++)
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
          for (j = 0; j < nonodes; j++) {
            bestyet = -10.0 * spp * chars;
            item = treenode[j];
            there = root;
            if (item != root) {
              re_move(item, &nufork, &root, recompute, treenode, &grbg, zeros);
              there = root;
              addpreorder(root, item, nufork);
              add(there, item, nufork, &root, recompute, treenode, &grbg, zeros);
            }
            if (progress) {
              putchar('.');
              fflush(stdout);
            }
          }
          if (progress)
            putchar('\n');
        } while (bestlike > gotlike);
      }
    }
    if (progress)
      putchar('\n');
    for (i = spp - 1; i >= 1; i--)
      re_move(treenode[i], &dummy, &root, recompute, treenode, &grbg, zeros);
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
      recompute = false;
      for (i = 0; i <= (nextree - 2); i++) {
        root = treenode[0];
        add(treenode[0], treenode[1], treenode[spp], &root, recompute,
              treenode, &grbg, zeros);
        for (j = 3; j <= spp; j++)
          add(treenode[bestrees[i].btree[j - 1] - 1], treenode[j - 1],
            treenode[spp + j - 2], &root, recompute, treenode, &grbg, zeros);
        reroot(treenode[outgrno - 1], root);
        postorder(root);
        evaluate(root);
        printree(root, 1.0);
        describe();
        for (j = 1; j < spp; j++)
          re_move(treenode[j], &dummy, &root, recompute, treenode, &grbg, zeros);
      }
    }
  } else {
    openfile(&intree, INTREE, "input tree file", "r", progname, intreename);
    fscanf(intree, "%ld%*[^\n]", &numtrees);
    getc(intree);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n\n\n");
    }
    fsteps = (long **)Malloc(maxuser*sizeof(long *));
    for (j = 1; j <= maxuser; j++)
      fsteps[j - 1] = (long *)Malloc(endsite*sizeof(long));
    names = (boolean *)Malloc(spp*sizeof(boolean));
    nodep = NULL;
    maxsteps = 0.0;
    which = 1;
    while (which <= numtrees) {
      firsttree = true;
      nextnode = 0;
      haslengths = true;
      treeread(intree, &root, treenode, &goteof, &firsttree,
                 nodep, intreename, "DnaComp", &nextnode, &haslengths,
                 &grbg, initdnacompnode);
      if (treeprint)
        fprintf(outfile, "\n\n");
      for (j = 0; j < spp; j++)
        names[j] = false;
      initboolnames(root, names);
      for (j = 0; j < spp; j++)
        in_tree[j] = names[j];
      j = 1;
      while (!in_tree[j - 1])
	j++;
      mincomp(j);
      if (outgropt)
	reroot(treenode[outgrno - 1], root);
      postorder(root);
      evaluate(root);
      printree(root, 1.0);
      describe();
      which++;
    }
    FClose(intree);
    putc('\n', outfile);
    if (numtrees > 1 && chars > 1 ) {
      fprintf(outfile, "Tree    Compatible  Difference  Its S.D.");
      fprintf(outfile, "   Significantly worse?\n\n");
      if (numtrees > maxuser)
        num = maxuser;
      else
        num = numtrees;
      for (which = 1; which <= num; which++) {
        fprintf(outfile, "%3ld%15.1f", which, nsteps[which - 1]);
        if (maxwhich == which)
          fprintf(outfile, "  <------ best\n");
        else {
          sumw = 0.0;
          sum = 0.0;
          sum2 = 0.0;
          for (j = 0; j < endsite; j++) {
            if (weight[j] > 0) {
              wt = weight[j];
              sumw += wt;
              sum += fsteps[maxwhich - 1][j] - fsteps[which - 1][j];
              TEMP2 = fsteps[which - 1][j] - fsteps[maxwhich - 1][j];
              sum2 += TEMP2 * TEMP2 / wt;
            }
          }
          sd = sqrt(sumw / (sumw - 1.0) * (sum2 - (sum /sumw) * (sum / sumw)));
          fprintf(outfile, "%9.1f%12.4f",
                  maxsteps - nsteps[which - 1], sd);
          if (sum > 1.95996 * sd)
            fprintf(outfile, "           Yes\n");
          else
            fprintf(outfile, "            No\n");
        }
      }
      fprintf(outfile, "\n\n");
    }
    for (j = 1; j <= maxuser; j++)
      free(fsteps[j - 1]);
    free(fsteps);
    free(names);
  }
  if (jumb == njumble) {
    if (progress) {
      printf("Output written to output file\n\n");
      if (trout)
        printf("Trees also written onto file\n\n");
    }
  }
}  /* maketree */


void freerest()
{
  if (!usertree) {
    freenode(&temp);
    freenode(&temp1);
  }
  freegrbg(&grbg);
  if (ancseq)
    freegarbage(&garbage);
  free(zeros);
  freenodes(nonodes, treenode);
}  /*  freerest */


int main(int argc, Char *argv[])
{  /* DNA compatibility by uphill search */
  /* reads in spp, chars, and the data. Then calls maketree to
    construct the tree */
#ifdef MAC
  argc = 1;		/* macsetup("Dnacomp","");	*/
  argv[0]="Dnacomp";
#endif
#ifdef WIN32
  phySetConsoleAttributes();
  phyClearScreen();
#endif
  openfile(&infile,INFILE,"input file", "r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file", "w",argv[0],outfilename);
  mulsets = false;
  garbage = NULL;
  grbg = NULL;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  datasets = 1;
  firstset = true;
  doinit();
  if (trout)
    openfile(&outtree,OUTTREE,"output tree file", "w",argv[0],outtreename);
  for (ith = 1; ith <= datasets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if (progress)
        printf("Data set # %ld:\n\n", ith);
    }
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
    freerest();
  }
  freetree(nonodes, treenode);
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
  exxit(0);
}  /* DNA compatibility by uphill search */


