

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


#include "parameters.h"
#include "mddata.h"
#include "water.h"
#include <math.h>

/* return the value of a with the same sign as b */
#define sign(a,b)  (b < 0 ) ? ( (a < 0) ? a : -a) : ( (a < 0) ? -a : a) 

CSHIFT(XA,XB,XMA,XMB,XL,BOXH,BOXL)
  /* compute some relevant distances between the two input molecules to
     this routine. If they are greater than the cutoff radius, compute
     these distances as if one of the particles were at its mirror image
     (periodic boundary conditions).
     Used by the intermolecular interactions routines */
  
  double XA[], XB[], XL[];
  double BOXH, BOXL, XMA, XMB;
  
{
    int I;
    
    /* Compute some intermolecular distances */
    
    XL[0] = XMA-XMB;
    XL[1] = XMA-XB[0];
    XL[2] = XMA-XB[2];
    XL[3] = XA[0]-XMB;
    XL[4] = XA[2]-XMB;
    XL[5] = XA[0]-XB[0];
    XL[6] = XA[0]-XB[2];
    XL[7] = XA[2]-XB[0];
    XL[8] = XA[2]-XB[2];
    XL[9] = XA[1]-XB[1];
    XL[10] = XA[1]-XB[0];
    XL[11] = XA[1]-XB[2];
    XL[12] = XA[0]-XB[1];
    XL[13] = XA[2]-XB[1];
    
    /* go through all 14 distances computed */
    for (I = 0; I <  14; I++) {
        /* if the value is greater than the cutoff radius */
        if (fabs(XL[I]) > BOXH) {
            XL[I]  =  XL[I] -(sign(BOXL,XL[I]));
        }
    } /* for */
    
} /* end of subroutine CSHIFT */


/* Generated from ../Source/cshift.C */
