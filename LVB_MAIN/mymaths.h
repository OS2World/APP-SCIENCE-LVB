/* LVB (c) Copyright 2003-2006 by Daniel Barker.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** mymaths.h - interface for mymaths.c ********** */

/* $Id: mymaths.h,v 1.9 2006/02/06 19:55:46 db60 Exp $ */

#ifndef LVB_MYMATHS_H
#define LVB_MYMATHS_H

double exp_wrapper(double) /*@globals errno@*/ /*@modifies nothing@*/ ;
double log_wrapper(double) /*@globals errno@*/ /*@modifies nothing@*/ ;
double pow_wrapper(double, double) /*@globals errno@*/ /*@modifies nothing@*/ ;

#endif /* LVB_MYMATHS_H */
