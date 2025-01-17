

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


#include "global.h"
#include "mdvar.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"

PREDIC(C,NOR1,ProcID)  /* predicts new values for displacement and its five
                          derivatives using Gear's sixth-order
                          predictor-corrector method */
  double C[];
  int NOR1;       /* NOR1 = NORDER + 1 = 7 (for a sixth-order method) */
  unsigned ProcID;
{
    /*   this routine calculates predicted F(X), F'(X), F''(X), ... */
    
    int JIZ;
    int  JI;
    int  L;
    int func, i, j, k, dir, atom;
    double S;
    struct link *curr_ptr;
    struct list_of_boxes *curr_box;
    
    curr_box = my_boxes[ProcID];
    
    while (curr_box) {
        
        i = curr_box->coord[XDIR];  /* X coordinate of box */
        j = curr_box->coord[YDIR];  /* Y coordinate of box */
        k = curr_box->coord[ZDIR];  /* Z coordinate of box */
        
        /* Loop through the current box's molecules */
        
        curr_ptr = BOX[i][j][k].list;
        
        while (curr_ptr) {
            
            JIZ = 2;
            
            /* loop over F(X), F'(X), F''(X), etc. */
            
            for (func = 0; func < NORDER; func++) {
                for ( dir = 0; dir < NDIR; dir++)
                    for ( atom = 0; atom < NATOM; atom++ ) {
                        JI = JIZ;
                        /* sum over Taylor Series */
                        S = 0.0;
                        for ( L = func; L < NORDER; L++) {
                            S += C[JI] * curr_ptr->mol.F[L+1][dir][atom];
                            JI++;
                        } /* for L */
                        curr_ptr->mol.F[func][dir][atom] += S;
                    } /* for atom */
                JIZ += NOR1;
            } /* for func */
            curr_ptr = curr_ptr->next_mol;
        } /* while curr_ptr */
        curr_box = curr_box->next_box;
    } /* while curr_box */
    
} /* end of subroutine PREDIC */

CORREC(PCC,NOR1,ProcID)
  double PCC[];       /* the predictor-corrector constants */
  int NOR1;           /* NORDER + 1 = 7 for a sixth-order method) */
  unsigned ProcID; 
  
  /* corrects the predicted values */
  
{
    
    /*   This routine calculates corrected F(X), F'(X), F"(X),
     *   from corrected F(X) = predicted F(X) + PCC(1)*(FR-SD)
     *   where SD is predicted accl. F"(X) and FR is computed 
     *   accl. (force/mass) at predicted position
     */
    
    double Y;
    int i, j, k, dir, atom, func;
    struct link *curr_ptr;
    box_list *curr_box;
    
    curr_box = my_boxes[ProcID];
    
    while (curr_box) {
        
        i = curr_box->coord[XDIR];  /* X coordinate of box */
        j = curr_box->coord[YDIR];  /* Y coordinate of box */
        k = curr_box->coord[ZDIR];  /* Z coordinate of box */
        
        /* Loop through the current box's molecules */
        
        curr_ptr = BOX[i][j][k].list;
        while (curr_ptr) {
            
            for (dir = 0; dir < NDIR; dir++) {
                for (atom = 0; atom < NATOM; atom++) {
                    Y = curr_ptr->mol.F[FORCES][dir][atom] - 
                        curr_ptr->mol.F[ACC][dir][atom];
                    for ( func = 0; func < NOR1; func++) 
                        curr_ptr->mol.F[func][dir][atom] += PCC[func] * Y;   
                } /* for atom */		
            } /* for dir */
            curr_ptr= curr_ptr->next_mol;
        } /* while curr_ptr */
        curr_box = curr_box->next_box;
    } /* while curr_box */
    
} /* end of subroutine CORREC */


/* Generated from ../Source/predcor.C */
