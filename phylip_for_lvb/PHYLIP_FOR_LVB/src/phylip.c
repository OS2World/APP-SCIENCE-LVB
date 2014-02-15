#if LVB
#include "../../LVB_MAIN/lvb.h"
#endif	/* LVB */
#include "phylip.h"

/* version 3.6. (c) Copyright 1993-1997 by Joseph Felsenstein.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   and Dan Fineman.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/* PHYLIP source code modified for use with LVB.
 * Modifications (c) Copyright 2003-2006 by Daniel Barker. */

FILE *infile, *outfile, *intree, *intree2, *outtree, *weightfile, *catfile;
long spp, words, bits;
boolean ibmpc, ansi, tranvsp;
naym *nayme;                     /* names of species */

#ifdef WIN32
#include <windows.h>
/* for console code (clear screen, text color settings) */
CONSOLE_SCREEN_BUFFER_INFO savecsbi;
HANDLE hConsoleOutput;

void phyClearScreen();
void phySaveConsoleAttributes();
void phySetConsoleAttributes();
void phyRestoreConsoleAttributes();
void phyFillScreenColor();
#endif


int eoff(FILE *f)
{
    register int ch;

    if (feof(f))
        return 1;
    if (f == stdin)
        return 0;
    ch = getc(f);
    if (ch == EOF)
        return 1;
    ungetc(ch, f);
    return 0;
}  /*eoff*/


int eoln(FILE *f)
{
    register int ch;

    ch = getc(f);
    if (ch == EOF)
        return 1;
    ungetc(ch, f);
    return (ch == '\n');
}  /*eoln*/


int filexists(char *filename)
{
  FILE *fp;
  fp =fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}  /*filexists*/


void get_command_name (char *vektor, char *progname_without_path)
{

  /* Puts the name of the program into progname from vektor without
     the whole path */
  char *last_slash;

  /* Point to the last slash... */
  last_slash = strrchr (vektor, '/');

  if (last_slash)
    /* If there was a last slash, point to the character after it */
    last_slash++;
  else
    /* If not, point to the vector */
    last_slash = vektor;
  
  strncpy (progname_without_path, last_slash, strlen (last_slash));
  progname_without_path [strlen(last_slash)] = '\0';
}  /*get_command_name*/


void getstryng(char *fname)
{ /* read in a file name from stdin and take off newline if any */

  fname = fgets(fname, 200, stdin);
  if (strchr(fname, '\n') != NULL)
    *strchr(fname, '\n') = '\0';
} /* getstryng */


void openfile(FILE **fp,char *filename,char *filedesc,
	char *mode,char *application,char *perm)
{
  /* open a file, testing whether it exists etc. */
  FILE *of;
  char file[100];
  char filemode[2];
  char input[100];
  char ch;
  char progname_without_path[100];

  get_command_name (application, progname_without_path);

  strcpy(file,filename);
  strcpy(filemode,mode);
  while (1){
    if (filemode[0] == 'w' && filexists(file)){
      printf("\n%s: the file \"%s\" that you wanted to\n",
          progname_without_path, file);
      printf("     use as %s already exists.\n", filedesc);
      printf("     Do you want to Replace it, Append to it,\n");
      printf("     write to a new File, or Quit?\n");
      do {
        printf("     (please type R, A, F, or Q) \n");
#ifdef WIN32
        phyFillScreenColor();
#endif
        fgets(input, sizeof(input), stdin);
        ch  = input[0];
        uppercase(&ch);
      } while (ch != 'A' && ch != 'R' && ch != 'F' && ch != 'Q');
      if (ch == 'Q')
        exxit(-1);
      if (ch == 'A') {
        strcpy(filemode,"a");
        continue;
      }
      else if (ch == 'F') {
       file[0] = '\0';
          while (file[0] =='\0') {
            printf("Please enter a new file name> ");
#ifdef WIN32
            phyFillScreenColor();
#endif
            getstryng(file);
          }
       strcpy(filemode,"w");
       continue;
        }
      }
    of = fopen(file,filemode);
    if (of)
      break;
    else {
      switch (filemode[0]){

      case 'r':
        printf("%s: can't find %s \"%s\"\n", progname_without_path,
            filedesc, file);
     file[0] = '\0';
        while (file[0] =='\0'){
          printf("Please enter a new file name> ");
#ifdef WIN32
          phyFillScreenColor();
#endif
          getstryng(file);}
        break;

      case 'w':
      case 'a':
        printf("%s: can't write %s file %s\n", progname_without_path,
            filedesc, file);
     file[0] = '\0';
        while (file[0] =='\0'){
          printf("Please enter a new file name> ");
#ifdef WIN32
          phyFillScreenColor();
#endif
          getstryng(file);}
        continue;
      default:
     printf("There is some error in the call of openfile. Unknown mode.\n");
     exxit(-1);
      }
    }
  }
  *fp=of;
  if (perm != NULL)
    strcpy(perm,file);
} /* openfile */


void cleerhome()
{
printf("%s", ((ibmpc || ansi) ? ("\033[2J\033[H") : "\n"));
} /* cleerhome */


double randum(longer seed)
{
  /* random number generator -- slow but machine independent */
  long i, j, k;
  long sum;
  longer mult, newseed;
  double x;

  mult[0] = 13;
  mult[1] = 24;
  mult[2] = 22;
  mult[3] = 6;
  for (i = 0; i <= 5; i++)
    newseed[i] = 0;
  for (i = 0; i <= 5; i++) {
    sum = newseed[i];
    k = i;
    if (i > 3)
      k = 3;
    for (j = 0; j <= k; j++)
      sum += mult[j] * seed[i - j];
    newseed[i] = sum;
    for (j = i; j <= 4; j++) {
      newseed[j + 1] += newseed[j] / 64;
      newseed[j] &= 63;
    }
  }
  memcpy(seed, newseed, sizeof(longer));
  seed[5] &= 3;
  x = 0.0;
  for (i = 0; i <= 5; i++)
    x = x / 64.0 + seed[i];
  x /= 4.0;
  return x;
}  /* randum */


void randumize(longer seed, long *enterorder)
{
  long i, j, k;

  for (i = 0; i < spp; i++) {
    j = (long)(randum(seed) * (i+1));
    k = enterorder[j];
    enterorder[j] = enterorder[i];
    enterorder[i] = k;
  }
} /* randumize */


long readlong(char *prompt)
{
long res;
char string[100];
do {
  printf("%s",prompt);
#ifdef WIN32
    phyFillScreenColor();
#endif
  getstryng(string);
  if (sscanf(string,"%ld",&res) == 1)
    break;
 } while (1);
return res;
}  /*readlong*/


void uppercase(Char *ch)
{
  /* convert ch to upper case */
  *ch = (islower (*ch) ? toupper(*ch) : (*ch));
}  /* uppercase */


void initseed(long *inseed, long *inseed0, longer seed)
{
  /* initialize random number seed */
  long i;

  printf("Random number seed (must be odd)?\n");
#ifdef WIN32
  phyFillScreenColor();
#endif
  scanf("%ld%*[^\n]", inseed);
  *inseed0 = *inseed;
  for (i = 0; i <= 5; i++)
    seed[i] = 0;
  i = 0;
  do {
    seed[i] = *inseed & 63;
    *inseed /= 64;
    i++;
  } while (*inseed != 0);
  getchar();
}  /*initseed*/


void initjumble(long *inseed, long *inseed0, longer seed, long *njumble)
{
  /* handle jumble option */
  initseed(inseed, inseed0, seed);
  printf("Number of times to jumble?\n");
#ifdef WIN32
  phyFillScreenColor();
#endif
  scanf("%ld%*[^\n]", njumble);
  getchar();
}  /*initjumble*/


void initoutgroup(long *outgrno, long spp)
{
  /* handle outgroup option */
  boolean done;

  done = true;
  do {
    printf("Type number of the outgroup:\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%ld%*[^\n]", outgrno);
    getchar();
    done = (*outgrno >= 1 && *outgrno <= spp);
    if (!done) {
      printf("BAD OUTGROUP NUMBER: %4ld\n", *outgrno);
      printf("  Must be in range 1 -%2ld\n", spp);
    }
  } while (done != true);
}  /*initoutgroup*/


void initthreshold(double *threshold)
{
  /* handle threshold option */
  boolean done;

  done = false;
  do {
    printf("What will be the threshold value?\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%lf%*[^\n]", threshold);
    getchar();
    done = (*threshold >= 1.0);
    if (!done)
      printf("BAD THRESHOLD VALUE:  it must be greater than 1\n");
    else
      *threshold = (long)(*threshold * 10.0 + 0.5) / 10.0;
  } while (done != true);
}  /*initthreshold*/


void initcatn(long *categs)
{ /* initialize category number */

  do {
    printf("Number of categories (1-%d)?\n", maxcategs);
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%ld%*[^\n]", categs);
    getchar();
  } while (*categs > maxcategs || *categs < 1);
}  /*initcatn*/


void initcategs(long categs, double *rate)
{ /* initialize category rates */
  long i;
  char line[100];
  char rest[100];
  long scanned;
  boolean done;

  for (;;){
    printf("Rate for each category? (use a space to separate)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    getstryng(line);
    done = true;
    for (i = 0; i < categs; i++){
      scanned = sscanf(line,"%lf %[^\n]", &rate[i],rest);
      if ((scanned < 2 && i < (categs - 1)) ||
       (scanned < 1 && i == (categs - 1))){
     printf("Please enter exactly %ld values.\n",categs);
     done = false;
     break;
      }
      strcpy(line,rest);
    }
    if (done)
      break;
  }
}  /*initcategs*/


void initprobcat(long categs, double *probsum, double *probcat)
{
  long i;
  boolean done;
  char line[100];
  char rest[100];
  long scanned;

  do {
    printf("Probability for each category?");
    printf(" (use a space to separate)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    getstryng(line);
    done = true;
    for (i = 0; i < categs; i++){
      scanned = sscanf(line,"%lf %[^\n]",&probcat[i],rest);
      if ((scanned < 2 && i < (categs - 1)) ||
       (scanned < 1 && i == (categs - 1))){
     done = false;
     printf("Please enter exactly %ld values.\n",categs);
     break;}
      strcpy(line,rest);
    }
    if (!done)
      continue;
    *probsum = 0.0;
    for (i = 0; i < categs; i++)
      *probsum += probcat[i];
    if (fabs(1.0 - (*probsum)) > 0.001) {
      done = false;
      printf("Probabilities must add up to");
      printf(" 1.0, plus or minus 0.001.\n");
    }
  } while (!done);
}  /*initprobcat*/



void inithowmany(long *howmany, long howoften)
{/* initialize how many cycles */
  do { 
    printf("How many cycles of %4ld trees?\n", howoften);
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%ld%*[^\n]", howmany);
    getchar();
  } while (*howmany <= 0);
}  /*inithowmany*/



void inithowoften(long *howoften)
{ /* initialize how many trees per cycle */
  do {
    printf("How many trees per cycle?\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%ld%*[^\n]", howoften);
    getchar();
  } while (*howoften <= 0);
}  /*inithowoften*/


void initlambda(double *lambda)
{
  do {
    printf("Mean block length of sites having the same rate (greater than 1)?\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%lf%*[^\n]", lambda);
    getchar();
    } while (*lambda <= 1.0);
  *lambda = 1.0 / *lambda;
}  /*initlambda*/


void initfreqs(double *freqa, double *freqc, double *freqg, double *freqt)
{
  char input[100];
  long scanned;
  printf("Base frequencies for A, C, G, T/U (use blanks to separate)?\n");
  do {
#ifdef WIN32
    phyFillScreenColor();
#endif
    getstryng(input);
    scanned = sscanf(input,"%lf%lf%lf%lf%*[^\n]", freqa, freqc, freqg, freqt);
    if (scanned == 4)
      break;
    else
      printf("Please enter exactly 4 values.\n");
  } while (1);
}  /*initfreqs*/


void initratio(double *ttratio)
{
  do {
    printf("Transition/transversion ratio?\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%lf%*[^\n]", ttratio);
    getchar();
  } while (*ttratio < 0.0);
}  /*initratio*/


void initpower(double *power)
{
  printf("New power?\n");
#ifdef WIN32
  phyFillScreenColor();
#endif
  scanf("%lf%*[^\n]", power);
  getchar();
}  /*initpower*/


void initdatasets(long *datasets)
{
  /* handle multi-data set option */
  boolean done;

  done = false;
  do {
    printf("How many data sets?\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%ld%*[^\n]", datasets);
    getchar();
    done = (*datasets >= 1);
      if (!done)
     printf("BAD DATA SETS NUMBER:  it must be greater than 1\n");
  } while (!done);
} /* initdatasets */

void justweights(long *datasets)
{
  /* handle multi-data set option by weights */
  boolean done;

  done = false;
  do {
    printf("How many sets of weights?\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%ld%*[^\n]", datasets);
    getchar();
    done = (*datasets >= 1);
      if (!done)
     printf("BAD NUMBER:  it must be greater than 1\n");
  } while (!done);
} /* justweights */


void initterminal(boolean *ibmpc, boolean *ansi)
{
  /* handle terminal option */
  if (*ibmpc) {
    *ibmpc = false;
    *ansi = true;
  } else if (*ansi)
      *ansi = false;
    else
      *ibmpc = true;
}  /*initterminal*/


void initnumlines(long *screenlines)
{
  do {
    *screenlines = readlong("Number of lines on screen?\n");
  } while (*screenlines <= 12);
}  /*initnumlines*/


void initbestrees(bestelm *bestrees, long maxtrees, boolean glob)
{
  /* initializes either global or local field of each array in bestrees */
  long i;

  if (glob)
    for (i = 0; i < maxtrees; i++)
      bestrees[i].gloreange = false;
  else
    for (i = 0; i < maxtrees; i++)
      bestrees[i].locreange = false;
} /* initbestrees */


void newline(FILE *filename, long i, long j, long k)
{
  /* go to new line if i is a multiple of j, indent k spaces */
  long m;

  if ((i - 1) % j != 0 || i <= 1)
    return;
  putc('\n', filename);
  for (m = 1; m <= k; m++)
    putc(' ', filename);
}  /* newline */


void inputnumbers(long *spp, long *chars, long *nonodes, long n)
{
  /* input the numbers of species and of characters */

  fscanf(infile, "%ld%ld", spp, chars);
  fscanf(infile, "%*[^\n]");
  *nonodes = *spp * 2 - n;
}  /* inputnumbers */


void inputnumbers2(long *spp, long *nonodes, long n)
{
  /* read species number */

  fscanf(infile, "%ld", spp);
  fscanf(infile, "%*[^\n]");
  fprintf(outfile, "\n%4ld Populations\n", *spp);
  *nonodes = *spp * 2 - n;
}  /* inputnumbers2 */


void inputnumbers3(long *spp, long *chars)
{
  /* input the numbers of species and of characters */

  fscanf(infile, "%ld%ld", spp, chars);
}  /* inputnumbers3 */


void samenumsp(long *chars, long ith)
{
  /* check if spp is same as the first set in other data sets */
  long cursp, curchs;

  if (eoln(infile)) {
    fscanf(infile, "%*[^\n]");
    getc(infile);
  }
  fscanf(infile, "%ld%ld", &cursp, &curchs);
  if (cursp != spp) {
    printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n", ith);
    exxit(-1);
  }
  *chars = curchs;
} /* samenumsp */


void samenumsp2(long ith)
{
  /* check if spp is same as the first set in other data sets */
  long cursp;

  if (eoln(infile)) {
    fscanf(infile, "%*[^\n]");
    getc(infile);
  }
  fscanf(infile, "%ld", &cursp);
  if (cursp != spp) {
    printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n", ith);
    exxit(-1);
  }
} /* samenumsp2 */


void readoptions(long *extranum, Char *options)
{ /* read option characters from input file */
  Char ch;

  while (!(eoln(infile))) {
    ch = getc(infile);
    uppercase(&ch);
    if (strchr(options, ch) != NULL)
     (* extranum)++;
    else if (ch != ' ') {
      printf("BAD OPTION CHARACTER: %c\n", ch);
      exxit(-1);
    }
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
}  /* readoptions */


void matchoptions(Char *ch, Char *options)
{  /* match option characters to those in auxiliary options line */

  *ch = getc(infile);
  uppercase(ch);
  if (strchr(options, *ch) == NULL) {
    printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE");
    printf(" WHICH STARTS WITH %c\n", *ch);
    exxit(-1);
  }
}  /* matchoptions */


void inputweights(long chars, steptr weight, boolean *weights)
{
  /* input the character weights, 0-9 and A-Z for weights 0 - 35 */
  Char ch;
  long i;

  for (i = 0; i < chars; i++) {
    do {
      if (eoln(weightfile)) {
        fscanf(weightfile, "%*[^\n]");
        getc(weightfile);
      }
      ch = getc(weightfile);
      if (ch == '\n')
        ch = ' ';
    } while (ch == ' ');
    weight[i] = 1;
    if (isdigit(ch))
      weight[i] = ch - '0';
    else if (isalpha(ch)) {
      uppercase(&ch);
      weight[i] = ch - 'A' + 10;
    } else {
      printf("BAD WEIGHT CHARACTER: %c\n", ch);
      exxit(-1);
    }
  }
  fscanf(weightfile, "%*[^\n]");
  getc(weightfile);
  *weights = true;
}  /* inputweights */


void inputweights2(long a, long b, long *weightsum,
	steptr weight, boolean *weights, char *prog)
{
  /* input the character weights, 0 or 1 */
  Char ch;
  long i;

  *weightsum = 0;
  for (i = a; i < b; i++) {
    do {
      if (eoln(weightfile)) {
        fscanf(weightfile, "%*[^\n]");
        getc(weightfile);
      }
      ch = getc(weightfile);
    } while (ch == ' ');
    weight[i] = 1;
    if (ch != '0' && ch != '1')
      weight[i] = ch - '0';
    else {
      printf("BAD WEIGHT CHARACTER: %c -- ", ch);
      printf("WEIGHTS IN %s MUST BE 0 OR 1\n", prog);
      exxit(-1);
    }
    *weightsum += weight[i];
  }
  *weights = true;
  fscanf(weightfile, "%*[^\n]");
  getc(weightfile);
}  /* inputweights2 */


void printweights(FILE *filename, long inc, long chars,
	steptr weight, Char *letters)
{
  /* print out the weights of sites */
  long i, j;

  fprintf(filename, "\n    %s are weighted as follows:\n",letters);
  for (i = 0; i < chars; i++) {
    if (i % 60 == 0) {
      putc('\n', filename);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', filename);
    }
    fprintf(filename, "%ld", weight[i + inc]);
    if ((i+1) % 10 == 0 && (i+1) % 60 != 0)
      putc(' ', filename);
  }
  fprintf(filename, "\n\n");
}  /* printweights */


void inputcategs(long a, long b, steptr category, long categs, char *prog)
{
  /* input the categories, 1-9 */
  Char ch;
  long i;

  for (i = a; i < b; i++) {
    do {
      if (eoln(catfile)) {
        fscanf(catfile, "%*[^\n]");
        getc(catfile);
      }
      ch = getc(catfile);
    } while (ch == ' ');
    if ((ch >= '1') && (ch <= ('0'+categs)))
      category[i] = ch - '0';
    else {
 printf("BAD CATEGORY CHARACTER: %c -- CATEGORIES IN %s ARE CURRENTLY 1-%ld\n",
             ch, prog, categs);
      exxit(-1);
    }
  }
  fscanf(catfile, "%*[^\n]");
  getc(catfile);
}  /* inputcategs */


void printcategs(FILE *filename, long chars, steptr category, Char *letters)
{
  /* print out the sitewise categories */
  long i, j;

  fprintf(filename, "\n    %s are:\n",letters);
  for (i = 0; i < chars; i++) {
    if (i % 60 == 0) {
      putc('\n', filename);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', filename);
    }
    fprintf(filename, "%ld", category[i]);
    if ((i+1) % 10 == 0 && (i+1) % 60 != 0)
      putc(' ', filename);
  }
  fprintf(filename, "\n\n");
}  /* printcategs */


void inputfactors(long chars, Char *factor, boolean *factors)
{
  /* reads the factor symbols */
  long i;
  Char ch;

  for (i = 1; i < nmlngth; i++) {
    ch = getc(infile);
    if (ch == '\n')
      ch = ' ';
  }
  for (i = 0; i < (chars); i++) {
    if (eoln(infile)) {
      fscanf(infile, "%*[^\n]");
      getc(infile);
    }
    factor[i] = getc(infile);
    if (factor[i] == '\n')
      factor[i] = ' ';
  }
  fscanf(infile, "%*[^\n]");
  getc(infile);
  *factors = true;
}  /* inputfactors */


void printfactors(FILE *filename, long chars, Char *factor, Char *letters)
{
  /* print out list of factor symbols */
  long i;

  fprintf(filename, "Factors%s:\n\n", letters);
  for (i = 1; i <= nmlngth - 5; i++)
    putc(' ', filename);
  for (i = 1; i <= (chars); i++) {
    newline(filename, i, 55, nmlngth + 3);
    putc(factor[i - 1], filename);
    if (i % 5 == 0)
      putc(' ', filename);
  }
  putc('\n', filename);
}  /* printfactors */


void headings(long chars, char *letters1, char *letters2)
{
  long i, j;

  putc('\n', outfile);
  j = nmlngth + (chars + (chars - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 37)
    j = 37;
  fprintf(outfile, "Name");
  for (i = 1; i <= j; i++)
    putc(' ', outfile);
  fprintf(outfile, "%s\n", letters1);
  fprintf(outfile, "----");
  for (i = 1; i <= j; i++)
    putc(' ', outfile);
  fprintf(outfile, "%s\n\n", letters2);
}  /* headings */


void initname(long i)
{
  /* read in species name */
  long j;

  for (j = 0; j < nmlngth; j++) {
    if (eoff(infile) | eoln(infile)){
      printf("ERROR: end-of-line or end-of-file");
      printf(" in the middle of species name for species %ld\n", i);
      exxit(-1);
    }
    nayme[i][j] = getc(infile);
  }
} /* initname */


void findtree(boolean *found,long *pos,long nextree,long *place,bestelm *bestrees)
{
  /* finds tree given by ARRAY place in ARRAY bestrees by binary search */
  /* used by dnacomp, dnapars, dollop, mix, & protpars */
  long i, lower, upper;
  boolean below, done;

  below = false;
  lower = 1;
  upper = nextree - 1;
  (*found) = false;
  while (!(*found) && lower <= upper) {
    (*pos) = (lower + upper) / 2;
    i = 3;
    done = false;
    while (!done) {
      done = (i > spp);
      if (!done)
        done = (place[i - 1] != bestrees[(*pos) - 1].btree[i - 1]);
      if (!done)
        i++;
    }
    (*found) = (i > spp);
    below = (place[i - 1] <  bestrees[(*pos )- 1].btree[i - 1]);
    if (*found)
      break;
    if (below)
      upper = (*pos) - 1;
    else
      lower = (*pos) + 1;
  }
  if (!(*found) && !below)
    (*pos)++;
}  /* findtree */


void addtree(long pos,long *nextree,boolean collapse,long *place,bestelm *bestrees)
{
  /* puts tree from ARRAY place in its proper position in ARRAY bestrees */
  /* used by dnacomp, dnapars, dollop, mix, & protpars */
  long i;

  for (i = *nextree - 1; i >= pos; i--){
    memcpy(bestrees[i].btree, bestrees[i - 1].btree, spp * sizeof(long));
    bestrees[i].gloreange = bestrees[i - 1].gloreange;
    bestrees[i - 1].gloreange = false;
    bestrees[i].locreange = bestrees[i - 1].locreange;
    bestrees[i - 1].locreange = false;
    bestrees[i].collapse = bestrees[i - 1].collapse;
  }
  for (i = 0; i < spp; i++)
    bestrees[pos - 1].btree[i] = place[i];
    bestrees[pos - 1].collapse = collapse;
  (*nextree)++;
}  /* addtree */


long findunrearranged(bestelm *bestrees, long nextree, boolean glob)
{
  /* finds bestree with either global or local field false */
  long i;

  if (glob) {
    for (i = 0; i < nextree - 1; i++)
      if (!bestrees[i].gloreange)
        return i;
  } else {
    for (i = 0; i < nextree - 1; i++)
      if (!bestrees[i].locreange)
        return i;
  }
  return -1;
} /* findunrearranged */


boolean torearrange(bestelm *bestrees, long nextree)
{ /* sees if any best tree is yet to be rearranged */

  if (findunrearranged(bestrees, nextree, true) >= 0)
    return true;
  else if (findunrearranged(bestrees, nextree, false) >= 0)
    return true;
  else
    return false;
} /* torearrange */


void reducebestrees(bestelm *bestrees, long *nextree)
{
  /* finds best trees with collapsible branches and deletes them */
  long i, j;

  i = 0;
  j = *nextree - 2;
  do {
    while (!bestrees[i].collapse && i < *nextree - 1) i++;
    while (bestrees[j].collapse && j >= 0) j--;
    if (i < j) {
      memcpy(bestrees[i].btree, bestrees[j].btree, spp * sizeof(long));
      bestrees[i].gloreange = bestrees[j].gloreange;
      bestrees[i].locreange = bestrees[j].locreange;
      bestrees[i].collapse = false;
      bestrees[j].collapse = true;
    }
  } while (i < j);
  *nextree = i + 1;
} /* reducebestrees */


void shellsort(double *a, long *b, long n)
{
  /* Shell sort keeping a, b in same order */
  /* used by dnapenny, dolpenny, & penny */
  long gap, i, j, itemp;
  double rtemp;

  gap = n / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= n; i++) {
      j = i - gap;
      while (j > 0) {
     if (a[j - 1] > a[j + gap - 1]) {
       rtemp = a[j - 1];
       a[j - 1] = a[j + gap - 1];
       a[j + gap - 1] = rtemp;
       itemp = b[j - 1];
       b[j - 1] = b[j + gap - 1];
       b[j + gap - 1] = itemp;
     }
     j -= gap;
      }
    }
    gap /= 2;
  }
}  /* shellsort */


void getch(Char *c, long *parens, FILE *treefile)
{
  /* get next nonblank character */
  do {
    if (eoln(treefile)) {
      fscanf(treefile, "%*[^\n]");
      (*c) = getc(treefile);
    }
    (*c) = getc(treefile);

    if ((*c) == '\n' || (*c) == '\t')
      (*c) = ' ';
  } while (!((*c) != ' ' || eoff(treefile)));
  if ((*c) == '(')
    (*parens)++;
  if ((*c) == ')')
    (*parens)--;
}  /* getch */


void getch2(Char *c, long *parens)
{
  /* get next nonblank character */
  do {
    if (eoln(intree)) {
      fscanf(intree, "%*[^\n]");
      getc(intree);
    }
    *c = getc(intree);
    if (*c == '\n' || *c == '\t')
      *c = ' ';
  } while (!(*c != ' ' || eoff(intree)));
  if (*c == '(')
   (*parens)++;
  if (*c == ')')
    (*parens)--;
}  /* getch2 */


void findch(Char c, Char *ch, long which)
{
  /* scan forward until find character c */
  boolean done;
  long dummy_parens;
  done = false;
  while (!done) {
    if (c == ',') {
      if (*ch == '(' || *ch == ')' || *ch == ';') {
        printf("\nERROR IN USER TREE%3ld: UNMATCHED PARENTHESIS OR MISSING COMMA\n",  which);
     exxit(-1);
      } else if (*ch == ',')
        done = true;
    } else if (c == ')') {
      if (*ch == '(' || *ch == ',' || *ch == ';') {
        printf("\nERROR IN USER TREE%3ld: ",which);
     printf("UNMATCHED PARENTHESIS OR NON-BIFURCATED NODE\n");
     exxit(-1);
      } else {
        if (*ch == ')')
          done = true;
      }
    } else if (c == ';') {
      if (*ch != ';') {
        printf("\nERROR IN USER TREE%3ld: ",which);
     printf("UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
     exxit(-1);
      } else
        done = true;
    }
    if (*ch != ')' && done)
      continue;
   getch(ch, &dummy_parens, intree);
  }
}  /* findch */


void findch2(Char c, long *lparens, long *rparens, Char *ch)
{
  /* skip forward in user tree until find character c */
  boolean done;
  long dummy_parens;
  done = false;
  while (!done) {
    if (c == ',') {
      if (*ch == '(' || *ch == ')' || *ch == ':' || *ch == ';') {
        printf("\nERROR IN USER TREE: ");
     printf("UNMATCHED PARENTHESIS OR MISSING COMMA\n");
        printf(" OR NON-TRIFURCATED BASE\n");
     exxit(-1);
      } else if (*ch == ',')
        done = true;
    } else if (c == ')') {
      if (*ch == '(' || *ch == ',' || *ch == ':' || *ch == ';') {
        printf("\nERROR IN USER TREE: UNMATCHED PARENTHESIS OR NON-BIFURCATED NODE\n");
     exxit(-1);
      } else if (*ch == ')') {
        (*rparens)++;
        if ((*lparens) > 0 && (*lparens) == (*rparens)) {
          if ((*lparens) == spp - 2) {
           getch(ch, &dummy_parens, intree);
            if (*ch != ';') {
              printf( "\nERROR IN USER TREE: ");
           printf("UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
           exxit(-1);
            }
          }
        }
     done = true;
      }
    }
    if (*ch != ')' && done)
      continue;
    if (*ch == ')')
     getch(ch, &dummy_parens, intree);
  }
}  /* findch2 */


void findch3(Char c, Char *ch, long which, long parens)
{
  /* scan forward until find character c */
  boolean done;
  double trweight;

  done = false;
  while (!done) {
    if (c == ';') {
      if ((*ch) == '[') {
        if (!eoln(intree)) {
          fscanf(intree, "%lf", &trweight);
          getc(intree);
          getch(ch, &parens, intree);
        }
      }
      if (*ch != ';') {
     printf("UNMATCHED PARENTHESIS OR MISSING SEMICOLON\n");
     exxit(-1);
      } else
        done = true;
    }
    if (*ch != ')' && done)
      continue;
    getch2(ch, &parens);
  }
  if (parens != 0) {
    printf("\nERROR IN TREE FILE:  UNMATCHED PARENTHESES\n");
    exxit(-1);
  }
}  /* findch3 */


void processlength(double *valyew, double *divisor, Char *ch, 
	boolean *minusread, FILE *treefile, long *parens)
{
  long digit, ordzero;
  boolean pointread;

  ordzero = '0';
  pointread = false;
  *minusread = false;
  *valyew = 0.0;
  *divisor = 1.0;
  getch(ch, parens, treefile);
  digit = (long)(*ch - ordzero);
  while ( ((digit <= 9) && (digit >= 0)) || *ch == '.' || *ch == '-') {
    if (*ch == '.' )
      pointread = true;
    else if (*ch == '-' )
      *minusread = true;
    else {
      *valyew = *valyew * 10.0 + digit;
      if (pointread)
        *divisor *= 10.0;
    }
    getch(ch, parens, treefile);
    digit = (long)(*ch - ordzero);
  }
  if (*minusread)
    *valyew = -(*valyew);
}  /* processlength */


void writename(long start, long n, long *enterorder)
{
  long i, j;

  for (i = start; i < start+n; i++) {
    printf(" %3ld. ", i+1);
    for (j = 0; j < nmlngth; j++)
      putchar(nayme[enterorder[i] - 1][j]);
    putchar('\n');
    fflush(stdout);
  }
}  /* writename */


void memerror()
{
printf("Error allocating memory\n");
exxit(-1);
}  /* memerror */


void odd_malloc(long x)
{
    printf ("ERROR: a function asked for an inappropriate amount of memory ");
    printf ("-- %ld bytes\n", x);
#if LVB
    printf ("       This can mean one of three things:\n");
    printf ("       1.  There is a bug in the program;\n");
    printf ("       2.  The data matrix is not a DNA matrix in PHYLIP 3.6 format; or\n");
    printf ("       3.  The input file is corrupt -- i.e. not saved as text-only.\n");
#else
    printf ("       This can mean one of two things:\n");
    printf ("       1.  There is a bug in the program, or\n");
    printf ("       2.  The input file is corrupt -- i.e. not saved as text-only.\n");
    printf ("\n");
    printf ("If it seems to be a bug, please mail joe@genetics.washington.edu\n");
    printf ("with a description of the problem and the input data file.\n");
#endif /* LVB */

    exxit(-1);
}


MALLOCRETURN *mymalloc(long x)
{
  MALLOCRETURN *new_block;

  if ((x <= 0) ||
      (x > TOO_MUCH_MEMORY))
    odd_malloc(x);

  new_block = (MALLOCRETURN *)calloc(1,x);

  if (!new_block) {
    memerror();
    return (MALLOCRETURN *) new_block;
  } else
    return (MALLOCRETURN *) new_block;
}

/*xx new code after this point-------------------------------*/


void gnu(node **grbg, node **p)
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */

  if (*grbg != NULL) {
    *p = *grbg;
    *grbg = (*grbg)->next;
  } else
    *p = (node *)Malloc(sizeof(node));

  (*p)->back       = NULL;
  (*p)->next       = NULL;
  (*p)->tip        = false;
  (*p)->times_in_tree = 0.0;
  (*p)->r          = 0.0;
  (*p)->theta      = 0.0;
  (*p)->x          = NULL;
  (*p)->protx	   = NULL;	/* for the sake of proml     */
}  /* gnu */


void chuck(node **grbg, node *p)
{
  /* collect garbage on p -- put it on front of garbage list */
  p->back = NULL;
  p->next = *grbg;
  *grbg = p;
}  /* chuck */


void zeronumnuc(node *p, long endsite)
{
  long i,j;

  for (i = 0; i < endsite; i++)
    for (j = (long)A; j <= (long)O; j++)
      p->numnuc[i][j] = 0;
} /* zeronumnuc */


void zerodiscnumnuc(node *p, long endsite)
{
  long i,j;

  for (i = 0; i < endsite; i++)
    for (j = (long)zero; j <= (long)seven; j++)
      p->discnumnuc[i][j] = 0;
} /* zerodiscnumnuc */


void allocnontip(node *p, long *zeros, long endsite)
{ /* allocate an interior node */
  /* used by dnacomp, dnapars, & dnapenny */

  p->numsteps = (steptr)Malloc(endsite*sizeof(long));
  p->oldnumsteps = (steptr)Malloc(endsite*sizeof(long));
  p->base = (baseptr)Malloc(endsite*sizeof(long));
  p->oldbase = (baseptr)Malloc(endsite*sizeof(long));
  p->numnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  memcpy(p->base, zeros, endsite*sizeof(long));
  memcpy(p->numsteps, zeros, endsite*sizeof(long));
  memcpy(p->oldbase, zeros, endsite*sizeof(long));
  memcpy(p->oldnumsteps, zeros, endsite*sizeof(long));
  zeronumnuc(p, endsite);
}  /* allocnontip */


void allocdiscnontip(node *p, long *zeros, unsigned char *zeros2, long endsite)
{ /* allocate an interior node */
  /* used by pars */

  p->numsteps = (steptr)Malloc(endsite*sizeof(long));
  p->oldnumsteps = (steptr)Malloc(endsite*sizeof(long));
  p->discbase = (discbaseptr)Malloc(endsite*sizeof(unsigned char));
  p->olddiscbase = (discbaseptr)Malloc(endsite*sizeof(unsigned char));
  p->discnumnuc = (discnucarray *)Malloc(endsite*sizeof(discnucarray));
  memcpy(p->discbase, zeros2, endsite*sizeof(unsigned char));
  memcpy(p->numsteps, zeros, endsite*sizeof(long));
  memcpy(p->olddiscbase, zeros2, endsite*sizeof(unsigned char));
  memcpy(p->oldnumsteps, zeros, endsite*sizeof(long));
  zerodiscnumnuc(p, endsite);
}  /* allocdiscnontip */


void allocnode(node **anode, long *zeros, long endsite)
{ /* allocate a node */
  /* used by dnacomp, dnapars, & dnapenny */

  *anode = (node *)Malloc(sizeof(node));
  allocnontip(*anode, zeros, endsite);
}  /* allocnode */


void allocdiscnode(node **anode, long *zeros, unsigned char *zeros2, 
	long endsite)
{ /* allocate a node */
  /* used by pars */

  *anode = (node *)Malloc(sizeof(node));
  allocdiscnontip(*anode, zeros, zeros2, endsite);
}  /* allocdiscnontip */


void gnutreenode(node **grbg, node **p, long i, long endsite, long *zeros)
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */

  if (*grbg != NULL) {
    *p = *grbg;
    *grbg = (*grbg)->next;
    memcpy((*p)->numsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->oldnumsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->base, zeros, endsite*sizeof(long));
    memcpy((*p)->oldbase, zeros, endsite*sizeof(long));
    zeronumnuc(*p, endsite);
  } else
    allocnode(p, zeros, endsite);
  (*p)->back = NULL;
  (*p)->next = NULL;
  (*p)->tip = false;
  (*p)->visited = false;
  (*p)->index = i;
  (*p)->numdesc = 0;
  (*p)->sumsteps = 0.0;
}  /* gnutreenode */


void gnudisctreenode(node **grbg, node **p, long i,
	long endsite, long *zeros, unsigned char *zeros2)
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */

  if (*grbg != NULL) {
    *p = *grbg;
    *grbg = (*grbg)->next;
    memcpy((*p)->numsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->oldnumsteps, zeros, endsite*sizeof(long));
    memcpy((*p)->discbase, zeros2, endsite*sizeof(unsigned char));
    memcpy((*p)->olddiscbase, zeros2, endsite*sizeof(unsigned char));
    zerodiscnumnuc(*p, endsite);
  } else
    allocdiscnode(p, zeros, zeros2, endsite);
  (*p)->back = NULL;
  (*p)->next = NULL;
  (*p)->tip = false;
  (*p)->visited = false;
  (*p)->index = i;
  (*p)->numdesc = 0;
  (*p)->sumsteps = 0.0;
}  /* gnudisctreenode */


void chucktreenode(node **grbg, node *p)
{
  /* collect garbage on p -- put it on front of garbage list */

  p->back = NULL;
  p->next = *grbg;
  *grbg = p;
}  /* chucktreenode */


void setupnode(node *p, long i)
{
  /* initialization of node pointers, variables */

  p->next = NULL;
  p->back = NULL;
  p->times_in_tree = (double) i * 1.0;
  p->index = i;
  p->tip = false;
}  /* setupnode */


long count_sibs (node *p)
{
  /* Count the number of nodes in a ring, return the total number of */
  /* nodes excluding the one passed into the function (siblings)     */
  node *q;
  long return_int = 0;

  if (p->tip) {
   printf ("Error: the function count_sibs called on a tip.  This is a bug.\n");
    exxit (-1);
  }

  q = p->next;
  while (q != p) {
    if (q == NULL) {
      printf ("Error: a loop of nodes was not closed.\n");
      exxit (-1);
    } else {
      return_int++;
      q = q->next;
    }
  }

  return return_int;
}  /* count_sibs */


void inittrav (node *p)
{
  /* traverse to set pointers uninitialized on inserting */
  long i, num_sibs;
  node *sib_ptr;
  
  if (p == NULL)
    return;
  if (p->tip)
    return;
  num_sibs = count_sibs (p);
  sib_ptr  = p;
  for (i=0; i < num_sibs; i++) {
    sib_ptr              = sib_ptr->next;
    sib_ptr->initialized = false;
    inittrav(sib_ptr->back);
  }
} /* inittrav */


void commentskipper(FILE ***intree, long *bracket)
{
  char c;
  
  c = fgetc(**intree);
  
  while (c != ']') {
    
    if(feof(**intree)) {
      printf("ERROR:  UNMATCHED COMMENT BRACKETS.");
      exxit(-1);
    }

    if(c == '[') {
      (*bracket)++;
      commentskipper(intree, bracket);
    }
    c = fgetc(**intree);
  }
  (*bracket)--;
}  /* commentskipper */


long countcomma(FILE **treefile, long *comma)
{
  /* Modified by Dan F. 11/10/96 */ 

  /* The next line inserted so this function leaves the file pointing
     to where it found it, not just re-winding it. */
  long orig_position = ftell (*treefile);

  Char c;
  long  lparen = 0;
  long bracket = 0;
  (*comma) = 0;


  for (;;){
    c = fgetc(*treefile);
    if (feof(*treefile))
      break;
    if (c == ';')
      break;
    if (c == ',')
      (*comma)++;
    if (c == '(')
         lparen++;
    if (c == '[') {
      bracket++;
      commentskipper(&treefile, &bracket);
    }
  }

  /* Don't just rewind, */
  /* rewind (*treefile); */
  /* Re-set to where it pointed when the function was called */

  fseek (*treefile, orig_position, SEEK_SET);

  return lparen + (*comma);
}  /*countcomma*/
/* countcomma rewritten so it passes back both lparen+comma to allocate nodep
   and a pointer to the comma variable.  This allows the tree to know how many
   species exist, and the tips to be placed in the front of the nodep array */


long countsemic(FILE **treefile)
{

  /* Used in dnaml to determine the number of user trees.  Return
     either a: the number of semicolons in the file outside comments
     or b: the first integer in the file */

  /* Re-written by Dan F 10/22/96 */
  Char c;
  long return_val, semic = 0;
  long bracket = 0;
  
  /* Eat all whitespace */
  c = fgetc (*treefile);
  while ((c == ' ')  ||
      (c == '\t') ||
      (c == '\n')) {
    c = fgetc (*treefile);
  }

  /* Then figure out if the first non-white character is a digit; if
     so, return it */
  if (isdigit (c)) {
    return_val = atoi(&c); 
  } else {

    /* Loop past all characters, count the number of semicolons
       outside of comments */
    for (;;){
      c = fgetc(*treefile);
      if (feof(*treefile))
     break;
      if (c == ';')
     semic++;
      if (c == '[') {
     bracket++;
     commentskipper(&treefile, &bracket);
      }
    }
    return_val = semic;
  }

  rewind (*treefile);
  return return_val;

}  /* countsemic */

#if 0
/* commenting out current countsemic() in favor of version
that Dan has rewritten.  That version takes just one arg, intree.  See
comments in that version.
Keeping this version in comments in case it is needed.  This change applied
to sccs version 1.4 -> 1.5
*/

long countsemic(FILE **intree, char *filename, char *application)
{
  Char c;
  long semic = 0;
  int bracket;

  bracket = 0;
  
  for (;;){
    c = fgetc(*intree);
    if (feof(*intree))
      break;
    if (c == ';')
      semic++;
    if (c == '[') {
      bracket++;
      commentskipper(&intree, &bracket);
    }
    
  }
  close(*intree);
  openfile(intree, INTREE ,"input tree file", "r", application, NULL);

  return semic;
}  /*countsemic;   Used in dnaml to determine the number of user trees.*/
#endif


void hookup(node *p, node *q)
{
  /* hook together two nodes */

  p->back = q;
  q->back = p;
}  /* hookup */


void link_trees(long local_nextnum, long nodenum, long local_nodenum,
	pointarray nodep)
{
  if(local_nextnum == 0)
    hookup(nodep[nodenum],nodep[local_nodenum]);
  
  else if(local_nextnum == 1)
    hookup(nodep[nodenum], nodep[local_nodenum]->next);
  
  else if(local_nextnum == 2)
    hookup(nodep[nodenum],nodep[local_nodenum]->next->next);
  
  else
    printf("Error in Link_Trees()");
} /* link_trees() */


void allocate_nodep(pointarray *nodep, FILE **treefile, long  *precalc_tips)  
{
  /* pre-compute space and allocate memory for nodep */

  long numnodes = 0;      /* returns number commas & (    */
  long numcom = 0;        /* returns number commas */
  
  numnodes = countcomma(treefile, &numcom) + 1;
  *nodep      = (pointarray)Malloc(2*numnodes*sizeof(node *));

  (*precalc_tips) = numcom + 1;        /* this will be used in placing the
       					tip nodes in the front region of
       					nodep.  Used for species check?  */
} /* allocate_nodep       -plc */


void malloc_pheno (node *p, long no_species, long endsite, long rcategs)
{
  /* Allocate the phenotype arrays; used by dnaml */
  long i;

  p->x  = (phenotype)Malloc(endsite*sizeof(ratelike));
  for (i = 0; i < endsite; i++)
    p->x[i]  = (ratelike)Malloc(rcategs*sizeof(sitelike));
} /* malloc_pheno */


void malloc_ppheno (node *p, long no_species, long endsite, long rcategs)
{
  /* Allocate the phenotype arrays; used by proml */
  long i;

  p->protx  = (pphenotype)Malloc(endsite*sizeof(pratelike));
  for (i = 0; i < endsite; i++)
    p->protx[i]  = (pratelike)Malloc(rcategs*sizeof(psitelike));
} /* malloc_ppheno */


long take_name_from_tree (Char *ch, Char *str, FILE *treefile)
{
  /* This loop takes in the name from the tree.
     Return the length of the name string.  */

  long name_length = 0;

  do {
    if ((*ch) == '_')
      (*ch) = ' ';
    str[name_length++] = (*ch);
    if (eoln(treefile)) {
      fscanf(treefile, "%*[^\n]");
      getc(treefile);
    }
    (*ch) = getc(treefile);
    if (*ch == '\n')
      *ch = ' ';
  } while ((*ch) != ':' && (*ch) != ',' && (*ch) != ')' &&
        (*ch) != '[' && (*ch) != ';' && name_length <= MAXNCH);
  return name_length;
}  /* take_name_from_tree */

void match_names_to_data (Char *str, pointarray treenode, node **p, long spp)
{
  /* This loop matches names taken from treefile to indexed names in
     the data file */

  boolean found;
  long i, n;

  n = 1;  
  do {
    found = true;
    for (i = 0; i < nmlngth; i++) {
      found = (found && ((str[i] == nayme[n - 1][i]) ||
        (((nayme[n - 1][i] == '_') && (str[i] == ' ')) ||
        ((nayme[n - 1][i] == ' ') && (str[i] == '\0')))));
    }
    
    if (found)
      *p = treenode[n - 1];
    else
      n++;

  } while (!(n > spp || found));
  
  if (n > spp) {
    printf("ERROR: CANNOT FIND SPECIES: ");
    for (i = 0; (str[i] != '\0') && (i < MAXNCH); i++)
      putchar(str[i]);
    putchar('\n');
    exxit(-1);
  }
}  /* match_names_to_data */


void addelement(node **p, node *q, Char *ch, long *parens, FILE *treefile,
	pointarray treenode, boolean *goteof, boolean *first, pointarray nodep,
	long *nextnode, long *ntips, boolean *haslengths, node **grbg,
	initptr initnode)
{

  /* recursive procedure adds nodes to user-defined tree */
  node *pfirst;
  long i, len = 0, nodei = 0;
  boolean notlast = 0;
  Char str[MAXNCH];
  node *r;

  if ((*ch) == '(') {

    (*nextnode)++;
    nodei = *nextnode;
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, bottom, treenode, nodep, str, ch, treefile);
    pfirst      = (*p);
    notlast = true;
    while (notlast) {
      (*initnode)(&(*p)->next, grbg, q,
                   len, nodei, ntips, parens, nonbottom, treenode,
                   nodep, str, ch, treefile);

      r = (*p)->next;
      getch(ch, parens, treefile);
      
      addelement(&(*p)->next->back, (*p)->next, ch, parens, treefile,
        treenode, goteof, first, nodep, nextnode, ntips,
        haslengths, grbg, initnode);

      (*initnode)(&r, grbg, q, len, nodei, ntips,
                    parens, hslength, treenode, nodep, str, ch, treefile);

      pfirst->numdesc++;
      *p = r;

      if ((*ch) == ')') {
        notlast = false;
        do {
          getch(ch, parens, treefile);
        } while ((*ch) != ',' && (*ch) != ')' &&
           (*ch) != '[' && (*ch) != ';' && (*ch) != ':');
      }
    }
    
    (*p)->next = pfirst;
    (*p)       = pfirst;

  } else if ((*ch) != ')') {
    for (i = 0; i < MAXNCH; i++) 
      str[i] = '\0';

    len = take_name_from_tree (ch, str, treefile);

    if ((*ch) == ')')
      (*parens)--;
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, tip, treenode, nodep, str, ch, treefile);
  } else
    getch(ch, parens, treefile);
  if (q != NULL)
    hookup(q, (*p));
  (*initnode)(p, grbg, q, len, nodei, ntips, 
                parens, iter, treenode, nodep, str, ch, treefile);
  if ((*ch) == ':')
    (*initnode)(p, grbg, q, len, nodei, ntips, 
                  parens, length, treenode, nodep, str, ch, treefile);
  else if ((*ch) != ';' && (*ch) != '[')
    (*initnode)(p, grbg, q, len, nodei, ntips, 
                  parens, hsnolength, treenode, nodep, str, ch, treefile);
  if ((*ch) == '[')
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, treewt, treenode, nodep, str, ch, treefile);
  else if ((*ch) == ';')
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, unittrwt, treenode, nodep, str, ch, treefile);
}  /* addelement */


void treeread (FILE *treefile, node **root, pointarray treenode,
	boolean *goteof, boolean *first, pointarray nodep, char *trefilename,
	char *application, long *nextnode, boolean *haslengths, node **grbg,
	initptr initnode)
{
  /* read in user-defined tree and set it up */
  char  ch;
  long parens = 0;
  long ntips = 0;
  
  (*goteof) = false;
  (*nextnode) = spp;

  while (eoln(treefile) && !eoff(treefile)) {
    /* Eats all whitespace (newlines and/or spaces) at the start of
       the treefile */
    fscanf(treefile, "%*[^\n]");
    getc(treefile);
  }

  if (eoff(treefile)) {
    (*goteof) = true;
    return;
  } 

  getch(&ch, &parens, treefile);

  while (ch != '(') {
    /* Eat everything in the file (i.e. digits, tabs) until you
       encounter an open-paren */
    getch(&ch, &parens, treefile);
  }
  (*haslengths) = true; 
  addelement(root, NULL, &ch, &parens, treefile,
         treenode, goteof, first, nodep, nextnode, &ntips,
         haslengths, grbg, initnode);
  while (eoln(treefile) && !eoff(treefile)) {
   /* Eats whitespace (newlines and/or spaces) between or after trees */
    fscanf(treefile, "%*[^\n]");
    getc(treefile);
  }
  (*first) = false;
  if (parens != 0) {
    printf("\nERROR IN TREE FILE:  UNMATCHED PARENTHESES\n");
    exxit(-1);
  }
}  /* treeread */


void addelement2(node *q, Char *ch, long *parens, FILE *treefile,
	pointarray treenode, boolean lngths, double *trweight, boolean *goteof,
	long *nextnode, long *ntips, long no_species, boolean *haslengths)
{
  /* recursive procedure adds nodes to user-defined tree
     -- old-style bifurcating-only version */
  node *pfirst = NULL, *p;
  long i, len, current_loop_index;
  boolean notlast, minusread;
  Char str[MAXNCH];
  double valyew, divisor;

  if ((*ch) == '(') {

    current_loop_index = (*nextnode) + spp;
    (*nextnode)++;

    /* This is an assignment of an interior node */
    p = treenode[current_loop_index];
    pfirst = p;
    notlast = true;
    while (notlast) {
      /* This while loop goes through a circle (triad for
      bifurcations) of nodes */
      p = p->next;
      /* added to ensure that non base nodes in loops have indices */
      p->index = current_loop_index + 1;
      
      getch(ch, parens, treefile);
      
      addelement2(p, ch, parens, treefile, treenode, lngths, trweight,
        goteof, nextnode, ntips, no_species, haslengths);

      if ((*ch) == ')') {
        notlast = false;
        do {
          getch(ch, parens, treefile);
        } while ((*ch) != ',' && (*ch) != ')' &&
           (*ch) != '[' && (*ch) != ';' && (*ch) != ':');
      }
    }

  } else if ((*ch) != ')') {
    for (i = 0; i < MAXNCH; i++) 
      str[i] = '\0';
    len = take_name_from_tree (ch, str, treefile);
    match_names_to_data (str, treenode, &p, spp);
    pfirst = p;
    if ((*ch) == ')')
      (*parens)--;
    (*ntips)++;
    strncpy (p->nayme, str, len);
  } else
    getch(ch, parens, treefile);
  
  if ((*ch) == '[') {    /* getting tree weight from last comment field */
    if (!eoln(treefile)) {
      fscanf(treefile, "%lf", trweight);
      getch(ch, parens, treefile);
      if (*ch != ']') {
        printf("ERROR: MISSING RIGHT SQUARE BRACKET\n");
        exxit(-1);
      }
      else {
        getch(ch, parens, treefile);
        if (*ch != ';') {
          printf("ERROR: MISSING SEMICOLON AFTER SQUARE BRACKETS\n");
          exxit(-1);
        }
      }
    }
  }
  else if ((*ch) == ';') {
    (*trweight) = 1.0 ;
    if (!eoln(treefile))
      printf("WARNING: TREE WEIGHT SET TO 1.0\n");
  }
  else
    (*haslengths) = ((*haslengths) && q == NULL);
  
  if (q != NULL)
    hookup(q, pfirst);
  if (q != NULL) {
    if (q->branchnum < pfirst->branchnum)
      pfirst->branchnum = q->branchnum;
    else
      q->branchnum = pfirst->branchnum;
  }

  if ((*ch) == ':') {
    processlength(&valyew, &divisor, ch,
       &minusread, treefile, parens);
    if (!minusread)
      p->oldlen = valyew / divisor;
    else
      p->oldlen = 0.0;
    if (lngths) {
      p->v = valyew / divisor;
      p->back->v = p->v;
      p->iter = false;
      p->back->iter = false;
    }
  }
  
}  /* addelement2 */


void treeread2 (FILE *treefile, node **root, pointarray treenode,
	boolean lngths, double *trweight, boolean *goteof,
	char *trefilename, char *application, boolean *haslengths,
        long *no_species)
{
  /* read in user-defined tree and set it up
     -- old-style bifurcating-only version */
  char  ch;
  long parens = 0;
  long ntips = 0;
  long nextnode;
  
  (*goteof) = false;
  nextnode = 0;

  while (eoln(treefile) && !eoff(treefile)) {
    /* Eats all whitespace (newlines and/or spaces) at the start of
       the treefile */
    fscanf(treefile, "%*[^\n]");
    getc(treefile);
  }

  if (eoff(treefile)) {
    (*goteof) = true;
    return;
  } 

  getch(&ch, &parens, treefile);

  while (ch != '(') {
    /* Eat everything in the file (i.e. digits, tabs) until you
       encounter an open-paren */
    getch(&ch, &parens, treefile);
  }

  addelement2(NULL, &ch, &parens, treefile, treenode, lngths, trweight,
          goteof, &nextnode, &ntips, (*no_species), haslengths);
  (*root) = treenode[*no_species];

  while (eoln(treefile) && !eoff(treefile)) {
    /* Eats whitespace (newlines and/or spaces) between or after trees */
    fscanf(treefile, "%*[^\n]");
    getc(treefile);
  }

  (*root)->oldlen = 0.0;

  if (parens != 0) {
    printf("\nERROR IN TREE FILE:  UNMATCHED PARENTHESES\n");
    exxit(-1);
  }
}  /* treeread2 */


void free_x (phenotype x, long endsite)
{
  /* Free the phenotype arrays of a given node */
  long i;

  for (i=0 ; i < endsite; i++) {
    if (x[i]) {
      free (x[i]);
      x[i] = NULL;
    } else
      continue;
  }
  free (x);
}

void free_x_of_ring (node *p, long endsite)
{
  /* Free all phenotype arrays in a loop (i.e. triad) */
  long i = 0;
  node *q;

  if (p->tip || (p->back == NULL)) {
    free_x (p->x, endsite);
    p->x = NULL;

  } else {

    q = p->next;
    while (q != p) {
      
      free_x (q->x, endsite);
      q->x = NULL;
      q    = q->next;
      
      if (++i > spp) {
     /* Something's wrong here, we've gone through this loop more
        times than there are species */
     printf ("Returning as the loop was done too many times\n");
     exxit (-1);
      }
    }
  }
}  /* free_x_of_ring */

void exxit(int exitcode)
{
  if (exitcode == 0)
    exit (exitcode);
  else {
#if LVB
    crash("see above");	/* error message will already have been given */
#else
    puts ("Hit Return to close program.  You may have to hit Return twice.");
    getchar ();
    getchar ();
#ifdef WIN32
    phyRestoreConsoleAttributes();
#endif
    exit (exitcode);
#endif	/* LVB */
  }
}

#ifdef WIN32
void phySaveConsoleAttributes()
{
  GetConsoleScreenBufferInfo( hConsoleOutput, &savecsbi );
}

void phySetConsoleAttributes()
{
  hConsoleOutput = GetStdHandle(STD_OUTPUT_HANDLE);

  phySaveConsoleAttributes();

  SetConsoleTextAttribute(hConsoleOutput, 
    BACKGROUND_GREEN | BACKGROUND_BLUE | BACKGROUND_INTENSITY);
}

void phyRestoreConsoleAttributes()
{
  COORD coordScreen = { 0, 0 };
  DWORD cCharsWritten;
  DWORD dwConSize; 

  dwConSize = savecsbi.dwSize.X * savecsbi.dwSize.Y;

  SetConsoleTextAttribute(hConsoleOutput, savecsbi.wAttributes);

  FillConsoleOutputAttribute( hConsoleOutput, savecsbi.wAttributes,
         dwConSize, coordScreen, &cCharsWritten );
}

void phyFillScreenColor()
{
  COORD coordScreen = { 0, 0 };
  DWORD cCharsWritten;
  CONSOLE_SCREEN_BUFFER_INFO csbi; /* to get buffer info */ 
  DWORD dwConSize; 

  GetConsoleScreenBufferInfo( hConsoleOutput, &csbi );
  dwConSize = csbi.dwSize.X * csbi.dwSize.Y;

  FillConsoleOutputAttribute( hConsoleOutput, csbi.wAttributes,
         dwConSize, coordScreen, &cCharsWritten );
}


void phyClearScreen()
{
   COORD coordScreen = { 0, 0 };    /* here's where we'll home the
                                       cursor */ 
   BOOL bSuccess;
   DWORD cCharsWritten;
   CONSOLE_SCREEN_BUFFER_INFO csbi; /* to get buffer info */ 
   DWORD dwConSize;                 /* number of character cells in
                                       the current buffer */ 

   /* get the number of character cells in the current buffer */ 

   GetConsoleScreenBufferInfo( hConsoleOutput, &csbi );
   dwConSize = csbi.dwSize.X * csbi.dwSize.Y;

   /* fill the entire screen with blanks */ 

   FillConsoleOutputCharacter( hConsoleOutput, (TCHAR) ' ',
      dwConSize, coordScreen, &cCharsWritten );

   /* get the current text attribute */ 

   GetConsoleScreenBufferInfo( hConsoleOutput, &csbi );

   /* now set the buffer's attributes accordingly */ 

   FillConsoleOutputAttribute( hConsoleOutput, csbi.wAttributes,
      dwConSize, coordScreen, &cCharsWritten );

   /* put the cursor at (0, 0) */ 

   SetConsoleCursorPosition( hConsoleOutput, coordScreen );
   return;
}
#endif

