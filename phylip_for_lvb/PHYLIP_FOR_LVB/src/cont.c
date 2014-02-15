#include "phylip.h"
#include "cont.h"

/* version 3.6. (c) Copyright 1999-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

extern long spp;


void alloctree(pointarray *treenode, long nonodes, boolean usertree)
{
  /* allocate treenode dynamically */
  /* used in contml & contrast */
  long i, j;
  node *p, *q;

  *treenode = (pointarray)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < spp; i++)
    (*treenode)[i] = (node *)Malloc(sizeof(node));
  for (i = spp; i < nonodes; i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    (*treenode)[i] = p;
  }
} /* alloctree */


void setuptree(tree *a, long nonodes, boolean usertree)
{
  /* initialize a tree */
  /* used in contml & contrast */
  long i, j;
  node *p;

  for (i = 1; i <= spp; i++) {
    a->nodep[i - 1]->back = NULL;
    a->nodep[i - 1]->tip = (i <= spp);
    a->nodep[i - 1]->iter = true;
    a->nodep[i - 1]->index = i;
  }
  for (i = spp + 1; i <= nonodes; i++) {
    p = a->nodep[i - 1];
    for (j = 1; j <= 3; j++) {
      p->back = NULL;
      p->tip = false;
      p->iter = true;
      p->index = i;
      p = p->next;
    }
  }
  a->likelihood = -99999.0;
  a->start = a->nodep[0];
}  /* setuptree */


void allocview(tree *a, long nonodes, long totalleles)
{
  /* allocate view */
  /* used in contml */
  long i, j;
  node *p;

  for (i = 0; i < spp; i++)
    a->nodep[i]->view = (phenotype3)Malloc(totalleles*sizeof(double));
  for (i = spp; i < nonodes; i++) {
    p = a->nodep[i];
    for (j = 1; j <= 3; j++) {
      p->view = (phenotype3)Malloc(totalleles*sizeof(double));
      p = p->next;
    }
  }
}  /* allocview */


void freeview(tree *a, long nonodes)
{
  /* deallocate view */
  /* used in contml */
  long i, j;
  node *p;

  for (i = 0; i < spp; i++)
    free(a->nodep[i]->view);
  for (i = spp; i < nonodes; i++) {
    p = a->nodep[i];
    for (j = 1; j <= 3; j++) {
      free(p->view);
      p = p->next;
    }
  }
}  /* freeview */

