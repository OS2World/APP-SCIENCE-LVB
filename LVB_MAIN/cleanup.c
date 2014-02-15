/* LVB (c) Copyright 2003-2006 by Daniel Barker.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** cleanup.c - prepare for exit ********** */

#include "lvb.h"

static const char *rcsid = "$Id: cleanup.c,v 1.20 2006/02/06 19:55:46 db60 Exp $";

Lvb_bool cleanup(void)
/* prevent apparent memory leaks to help debugging, log end time; return
 * LVB_TRUE on write error to stdout, LVB_FALSE otherwise */
{
    time_t endtim;	/* time at end of run */
    Lvb_bool val;	/* return value */

    endtim = time(NULL);
    printf("\n");
    printf("Ending at: %s", ctime(&endtim));
    printf("\n");

    /* log file won't be used again */
    fflush(stdout);
    if (ferror(stdout) != 0)
	val = LVB_TRUE;
    else
        val = LVB_FALSE;

    return val;
} /* end cleanup() */
