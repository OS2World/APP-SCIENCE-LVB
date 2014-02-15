/* LVB (c) Copyright 2003-2006 by Daniel Barker.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** selfconf.h ********** */

/* $Id: selfconf.h,v 1.8 2006/02/06 19:55:46 db60 Exp $ */

#ifndef LVB_SELFCONF_H
#define LVB_SELFCONF_H

typedef struct
{
    double t0;
    double t1;
    long maxpropose;
    long maxaccept;
    long maxfail;
    long runs;
    double predicted_length;
    double predicted_trees;
} Anneal_params_struct;

Anneal_params_struct self_configure(long) /*@globals internalState@*/;

#endif /* LVB_SELFCONF_H */
