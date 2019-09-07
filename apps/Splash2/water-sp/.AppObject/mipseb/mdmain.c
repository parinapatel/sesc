

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

#include "parameters.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "mddata.h"
#include "fileio.h"
#include "split.h"
#include "global.h"

/************************************************************************/

double  MDMAIN(NFSV,NFRST,NSTEP,NRST,NPRINT,NSAVE,LKT,NORD1,ProcID)
  int NFSV,NFRST,NSTEP,NRST,NPRINT,NSAVE,LKT,NORD1;
  unsigned ProcID;
  
{
    
    double TVIR = 0.0;
    double TTMV = 0.0;
    double TKIN = 0.0;
    double XTT;
    int i,j,k;
    double POTA,POTR,POTRF;
    double XVIR,AVGT,TEN;
    struct link *curr_ptr;
    struct list_of_boxes *new_box, *curr_box;
    
    for (i=start_end[ProcID]->box[XDIR][FIRST]; i<=start_end[ProcID]->box[XDIR][LAST]; i++) {
        for (j=start_end[ProcID]->box[YDIR][FIRST]; j<=start_end[ProcID]->box[YDIR][LAST]; j++) {
            for (k=start_end[ProcID]->box[ZDIR][FIRST]; k<=start_end[ProcID]->box[ZDIR][LAST]; k++) {
                new_box = (box_list *) our_malloc(sizeof(box_list),__FILE__,__LINE__);;
                new_box->coord[XDIR] = i;
                new_box->coord[YDIR] = j;
                new_box->coord[ZDIR] = k;
                new_box->next_box = NULL;
                curr_box = my_boxes[ProcID];
                if (curr_box == NULL)
                    my_boxes[ProcID] = new_box;
                else {
                    while (curr_box->next_box != NULL)
                        curr_box = curr_box->next_box;
                    curr_box->next_box = new_box;
                } /* else */
            }
        }
    }
    
    /* calculate initial value for acceleration */
    
    INTRAF(&gl->VIR,ProcID);
    
    {
 pthread_mutex_lock(&((gl->start).bar_mutex));
 (gl->start).bar_teller++;
 if ((gl->start).bar_teller == (NumProcs)) {
     (gl->start).bar_teller = 0;
     pthread_cond_broadcast(&((gl->start).bar_cond));
 } else
     pthread_cond_wait(&((gl->start).bar_cond), &((gl->start).bar_mutex));
 pthread_mutex_unlock(&((gl->start).bar_mutex));
};
    
    INTERF(ACC,&gl->VIR,ProcID);
    
    {
 pthread_mutex_lock(&((gl->start).bar_mutex));
 (gl->start).bar_teller++;
 if ((gl->start).bar_teller == (NumProcs)) {
     (gl->start).bar_teller = 0;
     pthread_cond_broadcast(&((gl->start).bar_cond));
 } else
     pthread_cond_wait(&((gl->start).bar_cond), &((gl->start).bar_mutex));
 pthread_mutex_unlock(&((gl->start).bar_mutex));
};
    
    /* MOLECULAR DYNAMICS LOOP */
    
    for (i=1;i <= NSTEP; i++) {
        TTMV=TTMV+1.00;
        
        /* POSSIBLE ENHANCEMENT:  Here's where one start measurements to avoid 
           cold-start effects.  Recommended to do this at the beginning of the
           second timestep; i.e. if (i == 2).
           */
        
        /* initialize various shared sums */
        if (ProcID == 0) {
            int dir;
            if (i >= 2) {
                {long time(); (gl->trackstart) = time(0);};
            }                
            gl->VIR = 0.0;
            gl->POTA = 0.0;
            gl->POTR = 0.0;
            gl->POTRF = 0.0;
            for (dir = XDIR; dir <= ZDIR; dir++)
                gl->SUM[dir] = 0.0;
        }
        
        if ((ProcID == 0) && (i >= 2)) {
            {long time(); (gl->intrastart) = time(0);};
        }
        
        {
 pthread_mutex_lock(&((gl->start).bar_mutex));
 (gl->start).bar_teller++;
 if ((gl->start).bar_teller == (NumProcs)) {
     (gl->start).bar_teller = 0;
     pthread_cond_broadcast(&((gl->start).bar_cond));
 } else
     pthread_cond_wait(&((gl->start).bar_cond), &((gl->start).bar_mutex));
 pthread_mutex_unlock(&((gl->start).bar_mutex));
};

        PREDIC(TLC,NORD1,ProcID);
        INTRAF(&gl->VIR,ProcID);
        
        {
 pthread_mutex_lock(&((gl->start).bar_mutex));
 (gl->start).bar_teller++;
 if ((gl->start).bar_teller == (NumProcs)) {
     (gl->start).bar_teller = 0;
     pthread_cond_broadcast(&((gl->start).bar_cond));
 } else
     pthread_cond_wait(&((gl->start).bar_cond), &((gl->start).bar_mutex));
 pthread_mutex_unlock(&((gl->start).bar_mutex));
};
        
        if ((ProcID == 0) && (i >= 2)) {
            {long time(); (gl->intraend) = time(0);};
            gl->intratime += gl->intraend - gl->intrastart;
        }
        
        if ((ProcID == 0) && (i >= 2)) {
            {long time(); (gl->interstart) = time(0);};
        }
        
        INTERF(FORCES,&gl->VIR,ProcID); 
        
        if ((ProcID == 0) && (i >= 2)) {
            {long time(); (gl->interend) = time(0);};
            gl->intertime += gl->interend - gl->interstart;
        }
        
        CORREC(PCC,NORD1,ProcID);
        
        BNDRY(ProcID);
        
        KINETI(NMOL,gl->SUM,HMAS,OMAS,ProcID);
        
        {
 pthread_mutex_lock(&((gl->start).bar_mutex));
 (gl->start).bar_teller++;
 if ((gl->start).bar_teller == (NumProcs)) {
     (gl->start).bar_teller = 0;
     pthread_cond_broadcast(&((gl->start).bar_cond));
 } else
     pthread_cond_wait(&((gl->start).bar_cond), &((gl->start).bar_mutex));
 pthread_mutex_unlock(&((gl->start).bar_mutex));
};
        
        if ((ProcID == 0) && (i >= 2)) {
            {long time(); (gl->intraend) = time(0);};
            gl->intratime += gl->intraend - gl->interend;
        }
        
        TKIN=TKIN+gl->SUM[0]+gl->SUM[1]+gl->SUM[2];
        TVIR=TVIR-gl->VIR;
        
        /* CHECK if  PRINTING AND/OR SAVING IS TO BE DONE */
        
        if ( ((i % NPRINT) == 0) || ((NSAVE > 0) && ((i % NSAVE) == 0))) {
            
            /* if so, call poteng to compute potential energy.  Note
               that we are attributing all the time in poteng to intermolecular
               computation although some of it is intramolecular (see poteng.C) */
            
            if ((ProcID == 0) && (i >= 2)) {
                {long time(); (gl->interstart) = time(0);};
            }
            
            POTENG(&gl->POTA,&gl->POTR,&gl->POTRF,ProcID);
            
            {
 pthread_mutex_lock(&((gl->start).bar_mutex));
 (gl->start).bar_teller++;
 if ((gl->start).bar_teller == (NumProcs)) {
     (gl->start).bar_teller = 0;
     pthread_cond_broadcast(&((gl->start).bar_cond));
 } else
     pthread_cond_wait(&((gl->start).bar_cond), &((gl->start).bar_mutex));
 pthread_mutex_unlock(&((gl->start).bar_mutex));
};
            
            if ((ProcID == 0) && (i >= 2)) {
                {long time(); (gl->interend) = time(0);};
                gl->intertime += gl->interend - gl->interstart;
            }
            
            POTA=gl->POTA*FPOT;
            POTR=gl->POTR*FPOT;
            POTRF=gl->POTRF*FPOT;
            XVIR=TVIR*FPOT*0.50/TTMV;
            AVGT=TKIN*FKIN*TEMP*2.00/(3.00*TTMV);
            TEN=(gl->SUM[0]+gl->SUM[1]+gl->SUM[2])*FKIN;
            XTT=POTA+POTR+POTRF+TEN;
            
            /* if it is time to print output as well ... */      
            if ((i % NPRINT) == 0 && ProcID == 0) {
                {pthread_mutex_lock(&(gl->IOLock));};
                fprintf(six,"     %5d %14.5lf %12.5lf %12.5lf %12.5lf \n"
                        ,i,TEN,POTA,POTR,POTRF);
                fprintf(six," %16.3lf %16.5lf %16.5lf\n",XTT,AVGT,XVIR);
                fflush(six);
                {pthread_mutex_unlock(&(gl->IOLock));};
            }
            
        }
        
        {
 pthread_mutex_lock(&((gl->start).bar_mutex));
 (gl->start).bar_teller++;
 if ((gl->start).bar_teller == (NumProcs)) {
     (gl->start).bar_teller = 0;
     pthread_cond_broadcast(&((gl->start).bar_cond));
 } else
     pthread_cond_wait(&((gl->start).bar_cond), &((gl->start).bar_mutex));
 pthread_mutex_unlock(&((gl->start).bar_mutex));
};
        
        if ((ProcID == 0) && (i >= 2)) {
            {long time(); (gl->trackend) = time(0);};
            gl->tracktime += gl->trackend - gl->trackstart;
        }
        
    } /* for i */
    
    return(XTT);
    
} /* mdmain.c */



/* Generated from ../Source/mdmain.C */
