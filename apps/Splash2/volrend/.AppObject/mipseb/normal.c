

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

/******************************************************************************
*                                                                             *
*   normal.c:  Compute normal map from density map, assuming gradient         *
*              operator is 1x3 and shading transition region is 0.            *
*                                                                             *
******************************************************************************/

#include "incl.h"

/* The following declarations show the layout of the .norm file.             */
/* If changed, the version number must be incremented and code               */
/* written to handle loading of both old and current versions.               */

				/* Version for new .norm files:              */
#define	NORM_CUR_VERSION   1	/*   Initial release                         */
short norm_version;		/* Version of this .norm file                */

short norm_len[NM];		/* Size of this normal map                   */

long norm_length;		/* Total number of normals in map            */
				/*   (= NM * product of lens)                */
NORMAL *norm_address;		/* Pointer to normal map                     */

/* End of layout of .norm file.                                              */

float nmag_epsilon;
void Normal_Compute();



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


#include "anl.h"

Compute_Normal()
{
  int i;

  for (i=0; i<NM; i++) {
    norm_len[i] = map_len[i]-2*INSET;
  }

  norm_length = norm_len[X] * norm_len[Y] * norm_len[Z];

  Allocate_Normal(&norm_address, norm_length);

  printf("    Computing normal...\n");

  Global->Index = NODE0;

#ifndef SERIAL_PREPROC
  for (i=1; i<num_nodes; i++) {
pthread_mutex_lock(&__intern__);
assert(__threads__<__MAX_THREADS__);
pthread_create(&(__tid__[__threads__++]), NULL, (void*(*)(void *))(Normal_Compute), NULL);
pthread_mutex_unlock(&__intern__);
}
#endif

  Normal_Compute();
}


Allocate_Normal(address, length)
     NORMAL **address;
     long length;
{
  unsigned int i,j,size,type_per_page,count,block;
  unsigned int p,numbytes;

  printf("    Allocating normal map of %ld bytes...\n",
	 length*sizeof(NORMAL));

  *address = (NORMAL *)our_malloc(length*sizeof(NORMAL),__FILE__,__LINE__);;

  if (*address == NULL)
    Error("    No space available for map.\n");

/*  POSSIBLE ENHANCEMENT:  Here's where one might distribute the 
    normal map among physical memories if one wanted to.
*/

  for (i=0; i<length; i++) *(*address+i) = 0;

}


void Normal_Compute()
{
  int inx,iny,inz;	/* Voxel location in object space            */
  int outx,outy,outz;	/* Loop indices in image space               */
  NORMAL *local_norm_address;
  double inv_mag_shift,magnitude,norm_lshift,grd_x,grd_y,grd_z;
  int nmap_partition,zstart,zstop;
  int xnorm,ynorm,znorm,norm0;
  int num_xqueue,num_yqueue,num_zqueue,num_queue;
  int xstart,xstop,ystart,ystop;
  int my_node;

  {pthread_mutex_lock(&(Global->IndexLock));};
  my_node = Global->Index++;
  {pthread_mutex_unlock(&(Global->IndexLock));};
  my_node = my_node%num_nodes;

/*  POSSIBLE ENHANCEMENT:  Here's where one might bind the process to a 
    processor, if one wanted to. 
*/

  /* use these to save work in loop */
  norm0 = (int)(NORM_LSHIFT);
  norm_lshift = -NORM_LSHIFT;
  nmag_epsilon = magnitude_epsilon;
  
  num_xqueue = ROUNDUP((float)norm_len[X]/(float)voxel_section[X]);
  num_yqueue = ROUNDUP((float)norm_len[Y]/(float)voxel_section[Y]);
  num_zqueue = ROUNDUP((float)norm_len[Z]/(float)voxel_section[Z]);
  num_queue = num_xqueue * num_yqueue * num_zqueue;
  xstart = (my_node % voxel_section[X]) * num_xqueue;
  xstop = MIN(xstart+num_xqueue,norm_len[X]);
  ystart = ((my_node / voxel_section[X]) % voxel_section[Y]) * num_yqueue;
  ystop = MIN(ystart+num_yqueue,norm_len[Y]);
  zstart = (my_node / (voxel_section[X] * voxel_section[Y])) * num_zqueue;
  zstop = MIN(zstart+num_zqueue,norm_len[Z]);

#ifdef SERIAL_PREPROC
  zstart=0;
  zstop = norm_len[Z];
  ystart=0;
  ystop = norm_len[Y];
  xstart=0;
  xstop = norm_len[X];
#endif

  for (outz=zstart; outz<zstop; outz++) {
    for (outy=ystart; outy<ystop; outy++) {
      for (outx=xstart; outx<xstop; outx++) {
	
	inx = INSET + outx;
	iny = INSET + outy;
	inz = INSET + outz;
	
	/* Compute voxel gradient assuming gradient operator is 1x3     */
	grd_x = (double)((int)MAP(inz,iny,inx+1) - (int)MAP(inz,iny,inx-1));
	grd_y = (double)((int)MAP(inz,iny+1,inx) - (int)MAP(inz,iny-1,inx));
	grd_z = (double)((int)MAP(inz+1,iny,inx) - (int)MAP(inz-1,iny,inx));

	/* Compute (magnitude*grd_divisor)**2 of gradient               */
	/* Reduce (magnitude*grd_divisor)**2 to magnitude*grd_divisor   */
	/* Reduce magnitude*grd_divisor to magnitude                    */
	magnitude = grd_x*grd_x+grd_y*grd_y+grd_z*grd_z; 

	local_norm_address = NORM_ADDRESS(outz,outy,outx,X);
	if (magnitude > SMALL) {
	  inv_mag_shift = norm_lshift/sqrt(magnitude); 
	  if (grd_x*inv_mag_shift > 0.0) xnorm = 1;
	  else xnorm = 0;
	  ynorm = (int)(grd_y*inv_mag_shift);
	  znorm = (int)(grd_z*inv_mag_shift);
	  *local_norm_address = TADDR(znorm,ynorm,xnorm);
	}
	else {
	  *local_norm_address = TADDR(norm0,2,1);
	}
      }
    }
  }

#ifndef SERIAL_PREPROC
  {
 pthread_mutex_lock(&((Global->SlaveBarrier).bar_mutex));
 (Global->SlaveBarrier).bar_teller++;
 if ((Global->SlaveBarrier).bar_teller == (num_nodes)) {
     (Global->SlaveBarrier).bar_teller = 0;
     pthread_cond_broadcast(&((Global->SlaveBarrier).bar_cond));
 } else
     pthread_cond_wait(&((Global->SlaveBarrier).bar_cond), &((Global->SlaveBarrier).bar_mutex));
 pthread_mutex_unlock(&((Global->SlaveBarrier).bar_mutex));
};
#endif
}


Load_Normal(filename)
     char filename[];
{
  char local_filename[FILENAME_STRING_SIZE];
  int fd;

  strcpy(local_filename,filename);
  strcat(local_filename,".norm");
  fd = Open_File(local_filename); 
  
  Read_Shorts(fd,&norm_version, (long)sizeof(norm_version));
  if (norm_version != NORM_CUR_VERSION) 
    Error("    Can't load version %d file\n",norm_version);
  
  Read_Shorts(fd,norm_len,(long)sizeof(map_len));
  
  Read_Longs(fd,&norm_length,(long)sizeof(norm_length));
  
  Allocate_Normal(&norm_address,norm_length);
  
  printf("    Loading normal map from .norm file...\n");
  Read_Shorts(fd,norm_address,(long)(norm_length*sizeof(NORMAL)));
  fflush(stdout);
  Close_File(fd);
}


Store_Normal(filename)
char filename[];
{
  char local_filename[FILENAME_STRING_SIZE];
  int fd;

  strcpy(local_filename,filename);
  strcat(local_filename,".norm");
  fd = Create_File(local_filename); 

  norm_version = NORM_CUR_VERSION;
  strcpy(local_filename,filename);
  strcat(local_filename,".norm");
  fd = Create_File(local_filename);
  Write_Shorts(fd,&norm_version,(long)sizeof(norm_version));
  
  Write_Shorts(fd,norm_len,(long)sizeof(norm_len));
  Write_Longs(fd,&norm_length,(long)sizeof(norm_length));
  
  printf("    Storing normal map into .norm file...\n");
  Write_Shorts(fd,norm_address,(long)(norm_length*sizeof(NORMAL)));
  Close_File(fd);
}


Deallocate_Normal(address)
NORMAL **address;
{
  printf("    Deallocating normal map...\n");

/*  G_FREE(*address);  */

  *address = NULL;
}


/* Generated from ../Source/normal.C */
