/* LVB (c) Copyright 2003-2006 by Daniel Barker.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** myuni.h - header for RNG functions in myuni.c ********** */

/* $Id: myuni.h,v 1.10 2006/02/06 19:55:46 db60 Exp $ */

#include <float.h>
#include <limits.h>

/* set max. random number seed value suitable for rinit() */
#if 900000001L > INT_MAX
#error LVB WARNING: type int not suitable, try with a 32-bit or larger system
#else
#define MAX_SEED 900000000
#endif  /* if 900000001L > INT_MAX */

/* external uni functions */
double uni(void);
void rinit(int ijkl);
