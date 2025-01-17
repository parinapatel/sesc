

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


#include <stdio.h>
#include "global.h"
#include "mdvar.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"


BNDRY(ProcID)     /* this routine puts the molecules back inside the box if
                     they are out */
  unsigned ProcID;
{
    int i, j, k, dir; 
    int X_INDEX, Y_INDEX, Z_INDEX;
    struct link *curr_ptr, *last_ptr, *next_ptr, *temp_ptr;
    struct list_of_boxes *curr_box;
    double *extra_p;
    
    /* for each box */
    curr_box = my_boxes[ProcID];
    while (curr_box) {
        i = curr_box->coord[XDIR];  /* X coordinate of box */
        j = curr_box->coord[YDIR];  /* Y coordinate of box */
        k = curr_box->coord[ZDIR];  /* Z coordinate of box */
        
        last_ptr = NULL;
        curr_ptr = BOX[i][j][k].list;
        
        /* Go through molecules in current box */
        
        while (curr_ptr) {  
            next_ptr = curr_ptr->next_mol;
            
            /* for each direction */
            
            for ( dir = XDIR; dir <= ZDIR; dir++ ) {
                extra_p = curr_ptr->mol.F[DISP][dir];
                
                /* if the oxygen atom is out of the box */
                if (extra_p[O] > BOXL ) {
                    
                    /* move all three atoms back in the box */
                    extra_p[H1] -= BOXL;
                    extra_p[O]  -= BOXL;
                    extra_p[H2] -= BOXL;
                }
                else if (extra_p[O] < 0.00) {
                    extra_p[H1] += BOXL;
                    extra_p[O]  += BOXL;
                    extra_p[H2] += BOXL;
                } 
            } /* for dir */
            
            /* If O atom moves out of current box, put it in correct box */
            X_INDEX = (int) (curr_ptr->mol.F[DISP][XDIR][O] / BOX_LENGTH);
            Y_INDEX = (int) (curr_ptr->mol.F[DISP][YDIR][O] / BOX_LENGTH);
            Z_INDEX = (int) (curr_ptr->mol.F[DISP][ZDIR][O] / BOX_LENGTH);
            
            if ((X_INDEX != i) || 
                (Y_INDEX != j) ||
                (Z_INDEX != k)) {
                
                /* Remove link from BOX[i][j][k] */
                
                {pthread_mutex_lock(&(BOX[i][j][k].boxlock));};
                if (last_ptr != NULL) 
                    last_ptr->next_mol = curr_ptr->next_mol;
                else
                    BOX[i][j][k].list = curr_ptr->next_mol;
                {pthread_mutex_unlock(&(BOX[i][j][k].boxlock));};
                
                /* Add link to BOX[X_INDEX][Y_INDEX][Z_INDEX] */
                
                {pthread_mutex_lock(&(BOX[X_INDEX][Y_INDEX][Z_INDEX].boxlock));};
                temp_ptr = BOX[X_INDEX][Y_INDEX][Z_INDEX].list;
                BOX[X_INDEX][Y_INDEX][Z_INDEX].list = curr_ptr;
                curr_ptr->next_mol = temp_ptr;
                {pthread_mutex_unlock(&(BOX[X_INDEX][Y_INDEX][Z_INDEX].boxlock));};
                
            }
            else last_ptr = curr_ptr;
            curr_ptr = next_ptr;   /* Go to next molecule in current box */
        } /* while curr_ptr */
        curr_box = curr_box->next_box;
    } /* for curr_box */
    
} /* end of subroutine BNDRY */


/* Generated from ../Source/bndry.C */
