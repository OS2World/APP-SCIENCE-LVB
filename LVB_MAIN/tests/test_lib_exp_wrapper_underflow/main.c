/* LVB (c) Copyright 2003-2006 by Daniel Barker.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include <lvb.h>

/* check exp_wrapper() crashes verbosely on underflow */

int main(void)
{
    lvb_initialize();

    exp_wrapper(-DBL_MAX);	/* should cause underflow */

    /* if we get this far, exp_wrapper() has failed to crash so the
     * test has failed */
    printf("test failed\n");

    exit(EXIT_SUCCESS);	/* we want failure: so program success
                         * indicates test failed */
}
