

/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/



#include <pthread.h>
#include <stdlib.h>
#include <semaphore.h>
#include <assert.h>
#if !(defined PAGE_SIZE)
#define PAGE_SIZE 4096
#endif
#define __MAX_THREADS__ 256

extern pthread_t __tid__[__MAX_THREADS__];
extern unsigned __threads__;
extern pthread_mutex_t __intern__;
void *our_malloc(size_t size, char * file, unsigned line);


#include <math.h>
#include "matrix.h"

#define MISS_COST 4.74
#define ALPHA 100
#define BETA 10
#define TIME(ops,misses) (ops+(misses)*MISS_COST)/(1.0+MISS_COST/BS)

extern BMatrix LB;
extern int P;
extern int BS;
double *opStats = NULL;
double seq_time, seq_ops, seq_misses;

PDIV(src_col, src_nz, ops, misses, runtime)
double *ops, *misses, *runtime;
{
  int super_size, passes;
  double this_ops, this_misses;

  super_size = src_col*src_nz - src_col*(src_col-1)/2;
  this_ops = 1.0*src_col*(src_col+1)*(2*src_col+1)/6;
  this_ops += 1.0*src_col*src_col*(src_nz-src_col);
  *ops += this_ops;

  passes = (src_col+BS-1)/BS;
  this_misses = 2.0*passes/3*super_size;
  *misses += this_misses;

  *runtime += TIME(this_ops, this_misses);
}


PMOD(src_col, dest_col, dest_nz, ops, misses, runtime)
double *ops, *misses, *runtime;
{
  int update_size, passes_src, passes_dest;
  double this_ops, this_misses;

  update_size = dest_col*dest_nz - dest_col*(dest_col-1)/2;
  this_ops = 2.0*src_col*update_size;
  *ops += this_ops;

  passes_src = (src_col+BS-1)/BS;
  passes_dest =  (dest_col+BS-1)/BS;
  if (passes_src == 1 && passes_dest == 1)
    /* just miss on source and dest once */
    this_misses = 1.0*src_col*dest_nz + 1.0*update_size;
  else
    this_misses = 1.0*passes_dest*src_col*dest_nz + 1.0*passes_src*update_size;
  *misses += this_misses;

  *runtime += TIME(this_ops, this_misses);
}

PADD(cols, rows, misses, runtime)
double *misses, *runtime;
{
  double this_misses;

  this_misses = 2*(rows*cols-cols*(cols-1)/2);
  *misses += this_misses;

  *runtime += TIME(0, this_misses);
}


AssignBlocksNow(distribute)
{
  int i, j, which;

  if (P == 1) {
    for (j=0; j<LB.n; j+=LB.partition_size[j])
      if (!LB.domain[j])
	for (i=LB.col[j]; i<LB.col[j+1]; i++)
	  BLOCK(i)->owner = 0;
  } else {
    EmbedBlocks(P);
  }
}


EmbedBlocks(P)
{
  int j, block;
  extern int scatter_decomposition;

  /* number the block partitions */

  for (j=0; j<LB.n; j+=LB.partition_size[j])
    if (!LB.domain[j]) {
      for (block=LB.col[j]; block<LB.col[j+1]; block++) {
	BLOCK(block)->owner = EmbeddedOwner(block);
      }
    }

  scatter_decomposition = 1;
}


/* Generated from ../Source/assign.C */
