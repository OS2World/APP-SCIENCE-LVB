/* LVB (c) Copyright 2003-2006 by Daniel Barker.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/**********

=head1 NAME

treestack.c - tree stack functions

Version Tag $Id: treestack.c,v 1.14 2006/02/06 19:55:47 db60 Exp $

=head1 DESCRIPTION

Provides operations for accessing and maintaining a stack of tree
topologies.

=cut

**********/

#include "lvb.h"

static const char *rcsid
 = "$Id: treestack.c,v 1.14 2006/02/06 19:55:47 db60 Exp $";

static void upsize(Treestack *sp)
/* increase allocation for tree stack *sp */
{
    extern Dataptr matrix;	/* data matrix */
    long i;	/* loop counter */

    sp->size++;
 
    /* allocate for stack itself */
    if (sp->stack == NULL)	/* 1st call, stack does not exist */
    {
        sp->stack = alloc(sp->size * sizeof(Treestack_element),
          "initial best tree stack");
        sp->next = 0;
        lvb_assert(sp->size == 1);	/* was incremented above */
    }
    else
    {
        sp->stack = realloc(sp->stack, sp->size * sizeof(Treestack_element));
        if (sp->stack == NULL)
            crash("out of memory: cannot increase allocation for\n"
             "best tree stack to %ld elements", sp->size);
    }

    /* allocate space within stack */
    for (i = sp->next; i < sp->size; i++)
    {
	sp->stack[i].tree = treealloc(brcnt(matrix->n), matrix->m);
	sp->stack[i].root = -1;
    }
 
} /* end upsize() */

static void dopush(Treestack *sp, const Branch *const barray, const long root)
/* push tree in barray (of root root) on to stack *sp */
{
    lvb_assert(sp->next <= sp->size);
    if (sp->next == sp->size)
        upsize(sp);
    treecopy(sp->stack[sp->next].tree, barray);
    sp->stack[sp->next].root = root;
    sp->next++;
 
} /* end dopush() */

/**********

=head1 treestack_cnt - RETURN COUNT OF TREES ON STACK

=head2 SYNOPSIS

long treestack_cnt(Treestack s);

=head2 DESCRIPTION

Return the number of trees currently stored on the stack C<s>.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item s

The tree stack whose size we wish to know.

=back

=head2 RETURN

Returns number of trees on stack C<s>.

=cut

**********/

long treestack_cnt(Treestack s)
{
    return s.next;

} /* end treestack_cnt() */

/**********

=head1 treestack_new - RETURN A NEW TREE STACK

=head2 SYNOPSIS

    Treestack treestack_new(void);

=head2 DESCRIPTION

Returns a new tree stack.

=head2 PARAMETERS

None.

=head2 RETURN

Returns a new, empty tree stack.

=cut

**********/ 

Treestack treestack_new(void)
{
    Treestack s;	/* return value */

    s.size = 0;
    s.next = 0;
    s.stack = NULL;

    return s;

} /* end treestack_new() */

/**********

=head1 treestack_push - PUSH TREE ONTO TREE STACK

=head2 SYNOPSIS

    long treestack_push(Treestack *sp, const Branch *const barray,
    const long root);

=head2 DESCRIPTION

Push copy of a tree onto an existing tree stack. Will not push if its
topology is already present on the stack. The stack will increase its
own memory allocation if necessary.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item barray

Pointer to first element of array containing tree to be pushed.

=item root

Root branch number of tree to be pushed.

=back

=head3 INOUT

=over 4

=item sp

Pointer to tree stack that will receive a new copy of the tree.

=back

=head2 RETURN

Returns 1 if the tree was pushed, or 0 if not.

=cut

**********/

long treestack_push(Treestack *sp, const Branch *const barray, const long root)
{
    long i;			/* loop counter */
    Branch *stacktree = NULL;	/* current tree on stack */
    long stackroot;		/* root of current tree */

    /* return before push if not a new topology */
    /* check backwards as similar trees may be discovered together */
    for (i = sp->next - 1; i >= 0; i--)
    {
        stacktree = sp->stack[i].tree;
        stackroot = sp->stack[i].root;
	if (treecmp(stacktree, stackroot, barray, root) == 0)
            return 0;
    }

    /* topology is new so must be pushed */
    dopush(sp, barray, root);
    return 1;

} /* end treestack_push() */

/**********

=head1 treestack_pop - POP TREE OFF TREE STACK

=head2 SYNOPSIS

    long treestack_pop(Branch *barray, long *root, Treestack *sp);
    
=head2 DESCRIPTION

Pop a tree off a tree stack.

=head2 PARAMETERS

=head3 OUTPUT

=over 4

=item barray

Pointer to first element of array to contain the popped tree. There
must be sufficient space in this array prior to the call.

=item root

Pointer to scalar that will receive the index of the root in C<barray>.

=back

=head3 INOUT

=over 4

=item sp

Pointer to tree stack from which we will pop the tree.

=back

=head2 RETURN

Returns 1 if a tree was popped, or 0 if the stack was empty.

=cut

**********/

long treestack_pop(Branch *barray, long *root, Treestack *sp)
{
    long val;	/* return value */

    if (sp->next >= 1)
    {
        sp->next--;
        treecopy(barray, sp->stack[sp->next].tree);
        *root = sp->stack[sp->next].root;

        val = 1;
    }
    else
    {
        val = 0;
    }

    return val;

} /* end treestack_pop() */

/**********

=head1 treestack_print - PRINT TREE STACK

=head2 SYNOPSIS

    long treestack_print(Treestack *sp, FILE *const outfp);

=head2 DESCRIPTION

Print all trees on a tree stack. The stack itself is not altered.

=head2 PARAMETERS

=head3 IN

=over 4

=item sp

The stack to be printed.

=back

=head3 OUTPUT

=over4

=item outfp

Pointer to the file to which we wish to output trees.

=back

=head2 RETURN

Returns the number of trees printed.

=cut

**********/

long treestack_print(Treestack *sp, FILE *const outfp)
{
    extern Dataptr matrix;	/* data matrix */
    const long d_obj1 = 0L;	/* 1st obj. for output trees */
    long root;			/* root of current tree */
    long i;			/* loop counter */
    Branch *barray;		/* current unpacked tree */

    /* "local" dynamic heap memory */
    barray = treealloc(brcnt(matrix->n), matrix->m);

    for (i = 0; i < sp->next; i++)
    {
        treecopy(barray, sp->stack[i].tree);
	root = objreroot(barray, sp->stack[i].root, d_obj1);
	lvb_treeprint(outfp, barray, root);
    }
    if (fflush(outfp) != 0)
	crash("file write error when writing best trees");

    /* deallocate "local" dynamic heap memory */
    free(barray);

    return sp->next;

} /* end bstprint() */

/**********

=head1 treestack_dump - DUMP AND CLEAR TREE STACK

=head2 SYNOPSIS

    long treestack_dump(Treestack *sp, FILE *const outfp);

=head2 DESCRIPTION

Pop all trees off a stack and dump them.

=head2 PARAMETERS

=head3 OUTPUT

=over4

=item outfp

Pointer to the file to which we wish to output trees.

=back

=head3 INOUT

=over 4

=item sp

The stack to be emptied and dumped.

=back

=head2 RETURN

Returns the number of trees dumped.

=cut

**********/

long treestack_dump(Treestack *sp, FILE *const outfp)
/* pop all trees on stack *sp and dump them to file outfp;
 * first branch (number 0); return number of trees dumped */
{
    extern Dataptr matrix;	/* data matrix */
    long cnt = 0;		/* tree count */
    long root;			/* number of root branch */
    Branch *barray;		/* current unpacked tree */

    /* "local" dynamic heap memory */
    barray = treealloc(brcnt(matrix->n), matrix->m);

    while ((treestack_pop(barray, &root, sp)) != 0)
    {
	treedump(outfp, barray);
	cnt++;
    }
    if (fflush(outfp) != 0)
	crash("file write error when dumping best trees");

    /* free "local" dynamic heap memory */
    free(barray);

    return cnt;

} /* end bstdump() */

/**********

=head1 treestack_free - DEALLOCATE TREE STACK

=head2 SYNOPSIS

    void treestack_free(Treestack *sp);

=head2 DESCRIPTION 

Clear a tree stack and deallocate dynamically allocated heap
memory associated with it.

=head2 PARAMETERS

=head3 INOUT 

=over 4

=item sp

The stack to be emptied and deallocated.

=back

=head2 RETURN

None.

=cut    

**********/

void treestack_free(Treestack *sp)
/* free all memory in tree stack *sp */
{
    long i;	/* loop counter */

    for (i = 0; i < sp->size; i++)
    {
	free(sp->stack[i].tree);
        sp->stack[i].tree = NULL;
        sp->stack[i].root = -1;
    }
    free(sp->stack);
    sp->next = 0;
    sp->size = 0;
    sp->stack = NULL;
 
} /* end bstfree() */

/**********

=head1 treestack_clear - EMPTY TREE STACK

=head2 SYNOPSIS

    void treestack_clear(Treestack *sp);

=head2 DESCRIPTION

Empty a tree stack but do not deallocate memory associated with it.
This memory will be available for re-use when trees are pushed onto the
stack again.

=head2 PARAMETERS

=head3 INOUT

=over 4

=item sp

The stack to be emptied.

=back

=head2 RETURN

None.

=cut

**********/

void treestack_clear(Treestack *sp)
/* clear stack *sp; note its allocation is not changed */
{
    sp->next = 0;	/* clear stack */

} /* end treestack_clear() */

/**********

=head1 treestack_transfer - TRANSFER TREES BETWEEN TREE STACKS

=head2 SYNOPSIS

    void treestack_transfer(Treestack *destp, Treestack *sourcep);

=head2 DESCRIPTION

Transfer one tree stack in its entirety to another. Order of the transferred
trees is not preserved. Current contents and current order of the destination
stack are preserved. The source stack is emptied but not deallocated.

=head2 PARAMETERS

=head3 INOUT

=over 4

=item destp

Pointer to the stack to be added to.

=item sourcep

Pointer to the stack to be transferred to C<destp> and cleared.

=back

=head2 RETURN

None.

=cut

**********/

void treestack_transfer(Treestack *destp, Treestack *sourcep)
{
    extern Dataptr matrix;	/* data matrix */
    Branch *barray;		/* current tree, in transit */
    long root;			/* number of root branch */

    /* "local" dynamic heap memory */
    barray = treealloc(brcnt(matrix->n), matrix->m);

    while (treestack_pop(barray, &root, sourcep) == 1)
    {
        treestack_push(destp, barray, root);
    }

    /* free "local" dynamic heap memory */
    free(barray);
}
