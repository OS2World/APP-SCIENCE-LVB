/* LVB (c) Copyright 2003-2006 by Daniel Barker.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include <lvb.h>

/* test for scream() */

int main(void)
{
    lvb_initialize();

    scream("%d %d %d %d %d %d %d %s!", 1, 2, 3, 4, 5, 6, 7,
    "TEST%STRING%GOES%d%d%dHERE");
    printf("returned!\n");
    return 0;
}
