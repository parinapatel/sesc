

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
*   octree.c:  Perform hiearchical enumeration of dataset.                    *
*                                                                             *
******************************************************************************/

#include "incl.h"

#define WRITE_PYR(IBIT,ILEVEL,IZ,IY,IX)\
				(PYR_ADDRESS(ILEVEL,IZ,IY,IX),\
				 *pyr_address2|=(IBIT)<<pyr_offset2)

/* The following declarations show the layout of the .pyr file.              */
/* If changed, the version number must be incremented and code               */
/* written to handle loading of both old and current versions.               */

				/* Version for new .pyr files:               */
#define	PYR_CUR_VERSION	1	/*   Initial release                         */
short pyr_version;		/* Version of this .pyr file                 */

short pyr_levels;		/* Number of levels in this pyramid          */
				
short pyr_len[MAX_PYRLEVEL+1][NM];	/* Number of voxels on each level    */
short pyr_voxlen[MAX_PYRLEVEL+1][NM];	/* Size of voxels on each level      */

long pyr_length[MAX_PYRLEVEL+1];/* Total number of bytes on this level       */
				/*   (= (int)((product of lens+7)/8))        */
BYTE *pyr_address[MAX_PYRLEVEL+1];/* Pointer to binary pyramid               */
				/*   (only pyr_levels sets of lens, lengths, */
				/*    and 3-D arrays are written to file)    */

/* End of layout of .pyr file.                                               */

				/* Subscripted access to binary pyramid      */
				/*   (PYR_ADDRESS only computes temporaries, */
				/*    most others invoke PYR_ADDRESS, then:  */
				/*      PYR returns char with bit in low bit,*/
				/*      CLEAR_PYR clears bit to 0,           */
				/*      SET_PYR sets bit to 1,               */
				/*      WRITE_PYR ORs BOOLEAN IBIT into bit, */
				/*      OVERWRITE_PYR clears, then ORs IBIT, */
				/*    CURRENT_PYR returns at current address)*/
				/* Warning:  do not invoke any of these      */
				/*   macros more than once per statement,    */
				/*   or values of temporaries may conflict,  */
				/*   depending on optimization by compiler.  */
long pyr_offset1,		/* Bit offset of desired bit within pyramid  */
     pyr_offset2;		/* Bit offset of bit within byte             */
BYTE *pyr_address2;		/* Pointer to byte containing bit            */

void Compute_Base();
void Or_Neighbors_In_Base();



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

Compute_Octree()
{
  int level,max_len;
  int i;
  max_len = 0;
  for (i=0; i<NM; i++) {
    max_len = MAX(max_len,opc_len[i]);
  }
  pyr_levels = 1;
  while ((1<<(pyr_levels-1)) < max_len)
    pyr_levels++;
  printf("    Computing binary pyramid of %d levels...\n",
	 pyr_levels);
  
  for (i=0; i<NM; i++) {
    pyr_len[0][i] = opc_len[i];
    pyr_voxlen[0][i] = 1;
  }
  pyr_length[0] = (pyr_len[0][X]*pyr_len[0][Y]*pyr_len[0][Z]+7)/8;
  Allocate_Pyramid_Level(&pyr_address[0], pyr_length[0]);

  printf("    Computing octree base...\n");

  Global->Index = NODE0;

/*  POSSIBLE ENHANCEMENT:  If you did want to replicate the octree
on all processors, don't do this create.
*/

#ifndef SERIAL_PREPROC
  for (i=1; i<num_nodes; i++) {
pthread_mutex_lock(&__intern__);
assert(__threads__<__MAX_THREADS__);
pthread_create(&(__tid__[__threads__++]), NULL, (void*(*)(void *))(Compute_Base), NULL);
pthread_mutex_unlock(&__intern__);
}
#endif

  Compute_Base();

  printf("    Performing OR of eight neighbors in binary mask...\n");

  Global->Index = 0;


/*  POSSIBLE ENHANCEMENT:  If you did want to replicate the octree
on all processors, don't do this create.
*/

#ifndef SERIAL_PREPROC
  for (i=1; i<num_nodes; i++) {
pthread_mutex_lock(&__intern__);
assert(__threads__<__MAX_THREADS__);
pthread_create(&(__tid__[__threads__++]), NULL, (void*(*)(void *))(Or_Neighbors_In_Base), NULL);
pthread_mutex_unlock(&__intern__);
}
#endif

  Or_Neighbors_In_Base();

  for (level=1; level<pyr_levels; level++) {
    for (i=0; i<NM; i++) {
      if (pyr_len[level-1][i] > 1) {
	pyr_len[level][i] =
	  (pyr_len[level-1][i]+1)>>1;
	pyr_voxlen[level][i] =
	  pyr_voxlen[level-1][i]<<1;
      }
      else {
	pyr_len[level][i] =
	  pyr_len[level-1][i];
	pyr_voxlen[level][i] =
	  pyr_voxlen[level-1][i];
      }
    }
    pyr_length[level] = (pyr_len[level][X]*
			 pyr_len[level][Y]*pyr_len[level][Z]+7)/8;
    Allocate_Pyramid_Level(&pyr_address[level], pyr_length[level]);
    Compute_Pyramid_Level(level);
  }
}


void Compute_Base()
{
  int outx, outy, outz, i;	
  int pmap_partition,zstart,zstop;
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

  num_xqueue = ROUNDUP((float)pyr_len[0][X]/(float)voxel_section[X]);
  num_yqueue = ROUNDUP((float)pyr_len[0][Y]/(float)voxel_section[Y]);
  num_zqueue = ROUNDUP((float)pyr_len[0][Z]/(float)voxel_section[Z]);
  num_queue = num_xqueue * num_yqueue * num_zqueue;
  xstart = (my_node % voxel_section[X]) * num_xqueue;
  xstop = MIN(xstart+num_xqueue,pyr_len[0][X]);
  ystart = ((my_node / voxel_section[X]) % voxel_section[Y]) * num_yqueue;
  ystop = MIN(ystart+num_yqueue,pyr_len[0][Y]);
  zstart = (my_node / (voxel_section[X] * voxel_section[Y])) * num_zqueue;
  zstop = MIN(zstart+num_zqueue,pyr_len[0][Z]);

/*  POSSIBLE ENHANCEMENT:  If you did want to replicate the octree
on all processors, then execute what's in the SERIAL_PREPROC ifdef below.
*/

#ifdef SERIAL_PREPROC
  zstart = 0;
  zstop = pyr_len[0][Z];
  ystart = 0;
  ystop = pyr_len[0][Y];
  xstart = 0;
  xstop = pyr_len[0][X];
#endif

  for (outz=zstart; outz<zstop; outz++) {
    for (outy=ystart; outy<ystop; outy++) {
      for (outx=xstart; outx<xstop; outx++) {
	if (OPC(outz,outy,outx) == 0)     /* mask bit is one unless opacity */
	  WRITE_PYR(0,0,outz,outy,outx);  /* value is zero.                 */
	else
	  WRITE_PYR(1,0,outz,outy,outx);
      }
    }
  }

/*  POSSIBLE ENHANCEMENT:  If you did want to replicate the octree
on all processors, then this computation is serial and you don't need
this barrier either.
*/

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


void Or_Neighbors_In_Base()
{
  int outx,outy,outz;	/* Loop indices in image space               */
  int outx_plus_one,outy_plus_one,outz_plus_one;
  int i;
  BOOLEAN bit;
  int pmap_partition,zstart,zstop;
  int my_node;

  {pthread_mutex_lock(&(Global->IndexLock));};
  my_node = Global->Index++;
  {pthread_mutex_unlock(&(Global->IndexLock));};
  my_node = my_node%num_nodes;

/*  POSSIBLE ENHANCEMENT:  Here's where one might bind the process to a 
    processor, if one wanted to. 
*/
    
  /* assumed for now that z direction has enough parallelism */
  pmap_partition = ROUNDUP((double)pyr_len[0][Z]/(double)num_nodes);
  zstart = pmap_partition * my_node;
  zstop = MIN(zstart+pmap_partition,pyr_len[0][Z]);

/*  POSSIBLE ENHANCEMENT:  If you did want to replicate the octree
on all processors, then you should execute what's in the SERIAL_PREPROC
  ifdef below.
*/

#ifdef SERIAL_PREPROC
  zstart = 0;
  zstop = pyr_len[0][Z];
#endif

  for (outz=zstart; outz<zstop; outz++) {
    outz_plus_one = MIN(outz+1,pyr_len[0][Z]-1);
    for (outy=0; outy<pyr_len[0][Y]; outy++) {
      outy_plus_one = MIN(outy+1,pyr_len[0][Y]-1);
      for (outx=0; outx<pyr_len[0][X]; outx++) {
	outx_plus_one = MIN(outx+1,pyr_len[0][X]-1);
	
	bit = PYR(0,outz,outy,outx);
	bit |= PYR(0,outz,outy,outx_plus_one);
	bit |= PYR(0,outz,outy_plus_one,outx);
	bit |= PYR(0,outz,outy_plus_one,outx_plus_one);
	bit |= PYR(0,outz_plus_one,outy,outx);
	bit |= PYR(0,outz_plus_one,outy,outx_plus_one);
	bit |= PYR(0,outz_plus_one,outy_plus_one,outx);
	bit |= PYR(0,outz_plus_one,outy_plus_one,outx_plus_one);
	
	WRITE_PYR(bit,0,outz,outy,outx);
      }
    }
  }

/*  POSSIBLE ENHANCEMENT:  If you did want to replicate the octree
on all processors, then this computation is serial and you don't
need this barrier either.
*/

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


Allocate_Pyramid_Level(address, length)
     BYTE **address;
     long length;
{
  unsigned int i,j,size,type_per_page,count,block;
  unsigned int p,numbytes;

  printf("      Allocating pyramid level of %ld bytes...\n",
	 length*sizeof(BYTE));

/*  POSSIBLE ENHANCEMENT:  If you want to replicate the octree
on all processors, then replace the macro below with a regular malloc.
*/

  *address = (BYTE *)our_malloc(length*sizeof(BYTE),__FILE__,__LINE__);;

  if (*address == NULL)
    Error("    No space available for pyramid level.\n");

/*  POSSIBLE ENHANCEMENT:  Here's where one might distribute the 
    octree among physical memories if one wanted to.
*/

  for (i=0; i<length; i++) *(*address+i) = 0;

}


Compute_Pyramid_Level(level)
     int level;
{
  int outx,outy,outz;	/* Loop indices in image space               */
  int i;
  int inx,iny,inz;
  int inx_plus_one,iny_plus_one,inz_plus_one;
  BOOLEAN bit;
  
  printf("      Computing pyramid level %d from level %d...\n",
	 level,level-1);
  for (outz=0; outz<pyr_len[level][Z]; outz++) {
    inz = outz<<1;
    inz_plus_one = MIN(inz+1,pyr_len[level-1][Z]-1);
    for (outy=0; outy<pyr_len[level][Y]; outy++) {
      iny = outy<<1;
      iny_plus_one = MIN(iny+1,pyr_len[level-1][Y]-1);
      for (outx=0; outx<pyr_len[level][X]; outx++) {
	inx = outx<<1;
	inx_plus_one = MIN(inx+1,pyr_len[level-1][X]-1);
	
	bit = PYR(level-1,inz,iny,inx);
	bit |= PYR(level-1,inz,iny,inx_plus_one);
	bit |= PYR(level-1,inz,iny_plus_one,inx);
	bit |= PYR(level-1,inz,iny_plus_one,inx_plus_one);
	bit |= PYR(level-1,inz_plus_one,iny,inx);
	bit |= PYR(level-1,inz_plus_one,iny,inx_plus_one);
	bit |= PYR(level-1,inz_plus_one,iny_plus_one,inx);
	bit |= PYR(level-1,inz_plus_one,iny_plus_one,inx_plus_one);
	
	WRITE_PYR(bit,level,outz,outy,outx);
      }
    }
  }
}


Load_Octree(filename)
     char filename[];
{
  char local_filename[FILENAME_STRING_SIZE];
  int fd,level;

  strcpy(local_filename,filename);
  strcat(local_filename,".pyr");
  fd = Open_File(local_filename);
  
  Read_Shorts(fd,&pyr_version, (long)sizeof(pyr_version));
  if (pyr_version != PYR_CUR_VERSION) 
    Error("    Can't load version %d file\n",pyr_version);
  
  Read_Shorts(fd,&pyr_levels,(long)sizeof(pyr_levels));
  Read_Shorts(fd,pyr_len,(long)(pyr_levels*NM*sizeof(short)));
  Read_Shorts(fd,pyr_voxlen,(long)(pyr_levels*NM*sizeof(short)));
  Read_Longs(fd,pyr_length,(long)(pyr_levels*sizeof(long)));

  printf("    Loading binary pyramid of %d levels...\n",pyr_levels);
  for (level=0; level<pyr_levels; level++) {
    Allocate_Pyramid_Level(&pyr_address[level],pyr_length[level]);
    printf("      Loading pyramid level %d from .pyr file...\n",level);
    Read_Bytes(fd,pyr_address[level],(long)(pyr_length[level]*sizeof(BYTE)));
  }
  Close_File(fd);
}


Store_Octree(filename)
char filename[];
{
  char local_filename[FILENAME_STRING_SIZE];
  int fd,level;

  strcpy(local_filename,filename);
  strcat(local_filename,".pyr");
  fd = Create_File(local_filename);

  pyr_version = PYR_CUR_VERSION;
  Write_Shorts(fd,&pyr_version,(long)sizeof(pyr_version));

  Write_Shorts(fd,&pyr_levels,(long)sizeof(pyr_levels));
  Write_Shorts(fd,pyr_len,(long)(pyr_levels*NM*sizeof(short)));
  Write_Shorts(fd,pyr_voxlen,(long)(pyr_levels*NM*sizeof(short)));
  Write_Longs(fd,pyr_length,(long)(pyr_levels*sizeof(long)));

  printf("    Storing binary pyramid of %d levels...\n",pyr_levels);
  for (level=0; level<pyr_levels; level++) {
    printf("      Storing pyramid level %d into .pyr file...\n",level);

    Write_Bytes(fd,pyr_address[level],(long)(pyr_length[level]*sizeof(BYTE)));
  }
  Close_File(fd);
}






/* Generated from ../Source/octree.C */
