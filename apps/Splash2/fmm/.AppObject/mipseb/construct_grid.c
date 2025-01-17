

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

#include <stdio.h>
#include <values.h>
#include "defs.h"
#include "memory.h"
#include "particle.h"
#include "box.h"
#include "partition_grid.h"
#include "construct_grid.h"

#define MY_PARTICLES     (Local[my_id].Particles)
#define MY_NUM_PARTICLES (Local[my_id].Num_Particles)
#define MY_MAX_PARTICLES (Local[my_id].Max_Particles)

void DetermineGridSize(int my_id);
void DetermineLocalGridSize(int my_id);
void MergeLocalGridSize(int my_id);
void ConstructLocalGrid(int my_id);
box *InitGrid(int my_id);
void InsertParticlesInTree(int my_id, particle **p_list, int num_of_particles,
			   box *root);
box *FindHome(int my_id, particle *p, box *current_home);
box *FindInitialRoot(int my_id, particle *p, box *current_home);
box *CreateChild(int my_id, box *pb, int new_child_num);
void SubdivideBox(int my_id, box *b);
void MergeLocalGrid(int my_id);
void MLGHelper(int my_id, box *local_box, box *global_box, box *global_parent);
void MergeLocalParticles(int my_id, particle **p_array, int num_of_particles,
			 box *pb);
void SplitParticles(int my_id, particle **p_array, int length,
		    particle **p_dist, int num_p_dist[NUM_OFFSPRING], box *pb);
box *CreateLeaf(int my_id, box *pb, int new_child_num, particle **p_array,
		int length);
void InsertParticlesInLeaf(int my_id, particle **p_array, int length, box *b);
int InsertBoxInGrid(int my_id, box *b, box *pb);
int RemoveBoxFromGrid(int my_id, box *b, box *pb);
void InsertSubtreeInPartition(int my_id, box *b);
void CleanupGrid(int my_id);
void SetSiblings(int my_id, box *b);
void SetColleagues(int my_id, box *b);
void ConstructGridLists(int my_id, box *b);
void ConstructInteractionLists(int my_id, box *b);
void SetVList(int my_id, box *b);
void SetUList(int my_id, box *b);
void SetUListHelper(int my_id, box *b, box *pb);
int AncestorBox(int my_id, box *b, box *ancestor_box);
void SetWList(int my_id, box *b);
void InsertNonAdjChildren(int my_id, box *b, box *pb);


void
ConstructGrid (int my_id, time_info *local_time, int time_all)
{
   unsigned long init, start, finish;
   int i;
  
   if (time_all)
      {long time(); (init) = time(0);};
   DetermineGridSize(my_id);   /* Finds the four corners of the grid. */
   FreeBoxes(my_id);
   InitPartition(my_id);
   if (time_all)
      {long time(); (start) = time(0);};
   if (MY_NUM_PARTICLES > 0) {
      ConstructLocalGrid(my_id);  /* Each processor constructs their own tree
				     based on only their particles */
      MergeLocalGrid(my_id);   /* The processors combine their trees into one
				  global tree. This step contains
				  communication between processors. */
   }
   {
 pthread_mutex_lock(&((G_Memory->synch).bar_mutex));
 (G_Memory->synch).bar_teller++;
 if ((G_Memory->synch).bar_teller == (Number_Of_Processors)) {
     (G_Memory->synch).bar_teller = 0;
     pthread_cond_broadcast(&((G_Memory->synch).bar_cond));
 } else
     pthread_cond_wait(&((G_Memory->synch).bar_cond), &((G_Memory->synch).bar_mutex));
 pthread_mutex_unlock(&((G_Memory->synch).bar_mutex));
};
   CleanupGrid(my_id);
   if (time_all)
      {long time(); (finish) = time(0);};

   if (time_all) {
      local_time[MY_TIME_STEP].other_time = start - init;
      local_time[MY_TIME_STEP].construct_time = finish - start;
   }	
}


void
ConstructLists (int my_id, time_info *local_time, int time_all)
{
   unsigned long start, finish;
  
   if (time_all)
      {long time(); (start) = time(0);};
   PartitionIterate(my_id, ConstructGridLists, TOP);
   {
 pthread_mutex_lock(&((G_Memory->synch).bar_mutex));
 (G_Memory->synch).bar_teller++;
 if ((G_Memory->synch).bar_teller == (Number_Of_Processors)) {
     (G_Memory->synch).bar_teller = 0;
     pthread_cond_broadcast(&((G_Memory->synch).bar_cond));
 } else
     pthread_cond_wait(&((G_Memory->synch).bar_cond), &((G_Memory->synch).bar_mutex));
 pthread_mutex_unlock(&((G_Memory->synch).bar_mutex));
};
   PartitionIterate(my_id, ConstructInteractionLists, BOTTOM);
   if (time_all)
      {long time(); (finish) = time(0);};

   if (time_all) {
      local_time[MY_TIME_STEP].list_time = finish - start;
   }	
}


void
DestroyGrid (int my_id, time_info *local_time, int time_all)
{
   box *b_scan, *tb;
   particle *p;
   int i;
   int particle_cost;
   unsigned long start, finish;
  
   if (time_all)
      {long time(); (start) = time(0);};
   b_scan = Local[my_id].Childless_Partition;
   MY_NUM_PARTICLES = 0;
   while (b_scan != NULL) {
      tb = b_scan;
      b_scan = b_scan->next;
      particle_cost = tb->cost / tb->num_particles;
      for (i = 0; i < tb->num_particles; i++) {
	 if (MY_MAX_PARTICLES <= MY_NUM_PARTICLES) {
	    LockedPrint("ERROR (P%d) : Too many particles in local array\n",
			my_id);
	    exit(-1);
	 }
	 p = tb->particles[i];
	 p->cost = particle_cost;
	 MY_PARTICLES[MY_NUM_PARTICLES++] = p;
      }
   }
   if (my_id == 0)
      Grid = NULL;
   if (time_all) {
      {long time(); (finish) = time(0);};
      local_time[MY_TIME_STEP].other_time += finish - start;
   }	
}


/*
 *  PrintGrid (int my_id)
 *
 *  Args : none.
 *
 *  Returns : nothing.
 *
 *  Side Effects : Prints the entire box structure of the grid to stdout.
 *
 */
void
PrintGrid (int my_id)
{
   if (Grid != NULL) {
      if (my_id == 0) {
	 printf("Info for Adaptive Grid :\n");
	 printf("Boxes :\n\n");
      }
      fflush(stdout);
      {
 pthread_mutex_lock(&((G_Memory->synch).bar_mutex));
 (G_Memory->synch).bar_teller++;
 if ((G_Memory->synch).bar_teller == (Number_Of_Processors)) {
     (G_Memory->synch).bar_teller = 0;
     pthread_cond_broadcast(&((G_Memory->synch).bar_cond));
 } else
     pthread_cond_wait(&((G_Memory->synch).bar_cond), &((G_Memory->synch).bar_mutex));
 pthread_mutex_unlock(&((G_Memory->synch).bar_mutex));
};
      PartitionIterate(my_id, PrintBox, TOP);
      {
 pthread_mutex_lock(&((G_Memory->synch).bar_mutex));
 (G_Memory->synch).bar_teller++;
 if ((G_Memory->synch).bar_teller == (Number_Of_Processors)) {
     (G_Memory->synch).bar_teller = 0;
     pthread_cond_broadcast(&((G_Memory->synch).bar_cond));
 } else
     pthread_cond_wait(&((G_Memory->synch).bar_cond), &((G_Memory->synch).bar_mutex));
 pthread_mutex_unlock(&((G_Memory->synch).bar_mutex));
};
      if (my_id == 0) {
	 printf("\n");
      }
      fflush(stdout);
      {
 pthread_mutex_lock(&((G_Memory->synch).bar_mutex));
 (G_Memory->synch).bar_teller++;
 if ((G_Memory->synch).bar_teller == (Number_Of_Processors)) {
     (G_Memory->synch).bar_teller = 0;
     pthread_cond_broadcast(&((G_Memory->synch).bar_cond));
 } else
     pthread_cond_wait(&((G_Memory->synch).bar_cond), &((G_Memory->synch).bar_mutex));
 pthread_mutex_unlock(&((G_Memory->synch).bar_mutex));
};
   }
   else
      printf("Adaptive grid has not been initialized yet.\n");
}


void
DetermineGridSize (int my_id)
{
   DetermineLocalGridSize(my_id);	/* Processor looks at its own particles and
				   finds the x and y max and min */
   MergeLocalGridSize(my_id);	/* Processors shares info with others and each
				   one computes the global max and min */
}


/* This code is pretty straightforward. A processor scans through its list of
   particles and finds the x_max, x_min, y_max, y_min. The only interesting
   thing is that particles are looked at two at a time, and compared to each
   other before looking at the current best max and min. This speeds up the
   running time of the algorithm from 2n to 3/2n. */
void
DetermineLocalGridSize (int my_id)
{
   real x_pos1, x_pos2, y_pos1, y_pos2;
   real x_max_challenger, x_min_challenger, y_max_challenger, y_min_challenger;
   int i;

   Local[my_id].Local_X_Max = -MAX_REAL;
   Local[my_id].Local_X_Min = MAX_REAL;
   Local[my_id].Local_Y_Max = -MAX_REAL;
   Local[my_id].Local_Y_Min = MAX_REAL;
   for (i = 0; i < MY_NUM_PARTICLES - 1; i += 2) {
      x_pos1 = MY_PARTICLES[i]->pos.x;
      y_pos1 = MY_PARTICLES[i]->pos.y;
      x_pos2 = MY_PARTICLES[i + 1]->pos.x;
      y_pos2 = MY_PARTICLES[i + 1]->pos.y;
      if (x_pos1 > x_pos2) {
	 x_max_challenger = x_pos1;
	 x_min_challenger = x_pos2;
      }
      else {
	 x_max_challenger = x_pos2;
	 x_min_challenger = x_pos1;
      }      
      if (y_pos1 > y_pos2) {
	 y_max_challenger = y_pos1;
	 y_min_challenger = y_pos2;
      }
      else {
	 y_max_challenger = y_pos2;
	 y_min_challenger = y_pos1;
      }
      if (x_max_challenger > Local[my_id].Local_X_Max)
	 Local[my_id].Local_X_Max = x_max_challenger;
      if (x_min_challenger < Local[my_id].Local_X_Min)
	 Local[my_id].Local_X_Min = x_min_challenger;
      if (y_max_challenger > Local[my_id].Local_Y_Max)
	 Local[my_id].Local_Y_Max = y_max_challenger;
      if (y_min_challenger < Local[my_id].Local_Y_Min)
	 Local[my_id].Local_Y_Min = y_min_challenger;
   }
   if (i == (MY_NUM_PARTICLES - 1)) {
      x_max_challenger = MY_PARTICLES[i]->pos.x;
      x_min_challenger = MY_PARTICLES[i]->pos.x;
      y_max_challenger = MY_PARTICLES[i]->pos.y;
      y_min_challenger = MY_PARTICLES[i]->pos.y;
      if (x_max_challenger > Local[my_id].Local_X_Max)
	 Local[my_id].Local_X_Max = x_max_challenger;
      if (x_min_challenger < Local[my_id].Local_X_Min)
	 Local[my_id].Local_X_Min = x_min_challenger;
      if (y_max_challenger > Local[my_id].Local_Y_Max)
	 Local[my_id].Local_Y_Max = y_max_challenger;
      if (y_min_challenger < Local[my_id].Local_Y_Min)
	 Local[my_id].Local_Y_Min = y_min_challenger;
   }
}


/* Each processor writes its best to a global 
   array, then they read everyone else's and find the absolute best. */
void
MergeLocalGridSize (int my_id)
{
   real *my_f_array, *their_f_array;
   real x_max_challenger, x_min_challenger, y_max_challenger, y_min_challenger;
   int i;
  
   my_f_array = G_Memory->f_array[my_id];
   my_f_array[0] = Local[my_id].Local_X_Max;
   my_f_array[1] = Local[my_id].Local_X_Min;
   my_f_array[2] = Local[my_id].Local_Y_Max;
   my_f_array[3] = Local[my_id].Local_Y_Min;
   {
 pthread_mutex_lock(&((G_Memory->synch).bar_mutex));
 (G_Memory->synch).bar_teller++;
 if ((G_Memory->synch).bar_teller == (Number_Of_Processors)) {
     (G_Memory->synch).bar_teller = 0;
     pthread_cond_broadcast(&((G_Memory->synch).bar_cond));
 } else
     pthread_cond_wait(&((G_Memory->synch).bar_cond), &((G_Memory->synch).bar_mutex));
 pthread_mutex_unlock(&((G_Memory->synch).bar_mutex));
};

   for (i = 0; i < Number_Of_Processors; i++) {
      their_f_array = G_Memory->f_array[i];
      x_max_challenger = their_f_array[0];
      x_min_challenger = their_f_array[1];
      y_max_challenger = their_f_array[2];
      y_min_challenger = their_f_array[3];
      if (x_max_challenger > Local[my_id].Local_X_Max)
	 Local[my_id].Local_X_Max = x_max_challenger;
      if (x_min_challenger < Local[my_id].Local_X_Min)
	 Local[my_id].Local_X_Min = x_min_challenger;
      if (y_max_challenger > Local[my_id].Local_Y_Max)
	 Local[my_id].Local_Y_Max = y_max_challenger;
      if (y_min_challenger < Local[my_id].Local_Y_Min)
	 Local[my_id].Local_Y_Min = y_min_challenger;
   }
}


void
ConstructLocalGrid (int my_id)
{
   Local[my_id].Local_Grid = InitGrid(my_id); /* Create the root box */
   InsertParticlesInTree(my_id, MY_PARTICLES, MY_NUM_PARTICLES, 
			 Local[my_id].Local_Grid);
   /* Put all of your particles into your local tree */
}


box *
InitGrid (int my_id)
{
   real x_length, y_length;
   real grid_length, grid_x_center, grid_y_center;
   int exp;
   box *ret_box;
  
   frexp(Local[my_id].Local_X_Max, &exp);
   if (Local[my_id].Local_X_Max > 0)
      Local[my_id].Local_X_Max = ldexp(1.0, exp);
   else {
      if (Local[my_id].Local_X_Max < 0)
	 Local[my_id].Local_X_Max = -ldexp(1.0, exp - 1);
   }
   frexp(Local[my_id].Local_X_Min, &exp);
   if (Local[my_id].Local_X_Min < 0)
      Local[my_id].Local_X_Min = -ldexp(1.0, exp);
   else {
      if (Local[my_id].Local_X_Min > 0)
	 Local[my_id].Local_X_Min = ldexp(1.0, exp - 1);
   }
   frexp(Local[my_id].Local_Y_Max, &exp);
   if (Local[my_id].Local_Y_Max > 0)
      Local[my_id].Local_Y_Max = ldexp(1.0, exp);
   else {
      if (Local[my_id].Local_Y_Max < 0)
	 Local[my_id].Local_Y_Max = -ldexp(1.0, exp - 1);
   }
   frexp(Local[my_id].Local_Y_Min, &exp);
   if (Local[my_id].Local_Y_Min < 0)
      Local[my_id].Local_Y_Min = -ldexp(1.0, exp);
   else {
      if (Local[my_id].Local_Y_Min > 0)
	 Local[my_id].Local_Y_Min = ldexp(1.0, exp - 1);
   }
  
   x_length = Local[my_id].Local_X_Max - Local[my_id].Local_X_Min;
   y_length = Local[my_id].Local_Y_Max - Local[my_id].Local_Y_Min;
   if (x_length > y_length)
      grid_length = x_length;
   else
      grid_length = y_length;
   grid_x_center = (grid_length / (real) 2.0) + Local[my_id].Local_X_Min;
   grid_y_center = (grid_length / (real) 2.0) + Local[my_id].Local_Y_Min;

   ret_box = InitBox(my_id, grid_x_center, grid_y_center, grid_length, NULL);
   return ret_box;
}


/* All this procedure does is put the particles in p_list into the tree pointed
   to by root. For each particle, you find out which childless (ie body) box
   the particle is supposed to go into, and insert it. Because childless boxes
   in this code can have 40 particles in them, instead of one for B-H, I have
   to check that. B-H is the same as having a MAX_PARTICLES_PER_BOX of 1. When
   that value is exceeded, the childless box is subdivided and the particles
   are divided up amongst the children. Note that all the particles may be put
   into one child, in which case that child must be subdivided as well, and so
   on until there is more than one child. */
void
InsertParticlesInTree (int my_id, particle **p_list, int num_of_particles, box *root)
{
   particle *p;
   box *dest_box;
   int i, j;
  
   dest_box = root;
   for (i = 0; i < num_of_particles; i++) {
      p = p_list[i];
      dest_box = FindHome(my_id, p, dest_box);
      dest_box->particles[dest_box->num_particles++] = p;
      while (dest_box->num_particles > MAX_PARTICLES_PER_BOX) {
	 SubdivideBox(my_id, dest_box);
	 if (dest_box->num_children == 1) {
	    for (j = 0; dest_box->children[j] == NULL; j++)
	       ;
	    dest_box = dest_box->children[j];
	 }
      }
   }
}


/* This function compares the particles position to the center of the parent 
   (or cell) box and chooses the appropriate child to move down to, until 
   either a child box or a null pointer is reached. In the first case, the box 
   is returned, and in the second, a new child box is created for the parent,
   and that box is returned. */
box *
FindHome (int my_id, particle *p, box *current_home)
{
   box *pb;

   pb = FindInitialRoot(my_id, p, current_home);
   while (pb->type == PARENT) {
      if (p->pos.y > pb->y_center) {
	 if (p->pos.x > pb->x_center) {
	    if (pb->children[0] == NULL)
	       CreateChild(my_id, pb, 0);
	    pb = pb->children[0];
	 }
	 else {
	    if (pb->children[1] == NULL)
	       CreateChild(my_id, pb, 1);
	    pb = pb->children[1];
	 }
      }
      else {
	 if (p->pos.x > pb->x_center) {
	    if (pb->children[3] == NULL)
	       CreateChild(my_id, pb, 3);
	    pb = pb->children[3];
	 }
	 else {
	    if (pb->children[2] == NULL)
	       CreateChild(my_id, pb, 2);
	    pb = pb->children[2];
	 }
      }
   }
   return pb;
}


box *
FindInitialRoot (int my_id, particle *p, box *current_home)
{
   int found;
   real x_center_distance, y_center_distance;

   found = FALSE;
   while (found == FALSE) {
      x_center_distance = p->pos.x - current_home->x_center;
      y_center_distance = p->pos.y - current_home->y_center;
      if (x_center_distance < 0)
	 x_center_distance = -x_center_distance;
      if (y_center_distance < 0)
	 y_center_distance = -y_center_distance;
      if ((x_center_distance > (current_home->length / 2.0)) ||
	  (y_center_distance > (current_home->length / 2.0)))
	 current_home = current_home->parent;
      else
	 found = TRUE;
   }
   return current_home;
}    
  


/* Simply creates a new box and sets the parent and child pointers correctly.
   The child's spot in the parent's children array determines its position.
   So, 0 is upper right, 1 is upper left, 2 is bottom left, and 3 is bottom
   right. That's how I know how to set the location of the child box's center
   just by knowing the child number. */
box *
CreateChild (int my_id, box *pb, int new_child_num)
{
   real child_length, child_offset;
   box *ret_box;
  
   child_length = pb->length / (real) NUM_DIMENSIONS;
   child_offset = pb->length / (real) NUM_OFFSPRING;  
   if (new_child_num == 0) {
      pb->children[0] = InitBox(my_id, (pb->x_center + child_offset), 
				(pb->y_center + child_offset), child_length, 
				pb);
      pb->shadow[0] = pb->children[0];
   }
   if (new_child_num == 1) {
      pb->children[1] = InitBox(my_id, (pb->x_center - child_offset), 
				(pb->y_center + child_offset), child_length, 
				pb);
      pb->shadow[1] = pb->children[1];
   }
   if (new_child_num == 2) {
      pb->children[2] = InitBox(my_id, (pb->x_center - child_offset), 
				(pb->y_center - child_offset), child_length, 
				pb);
      pb->shadow[2] = pb->children[2];
   }
   if (new_child_num == 3) {
      pb->children[3] = InitBox(my_id, (pb->x_center + child_offset), 
				(pb->y_center - child_offset), child_length, 
				pb);
      pb->shadow[3] = pb->children[3];
   }
   pb->children[new_child_num]->child_num = new_child_num;
   ret_box = pb->children[new_child_num];
   pb->num_children += 1;
   return ret_box;
}


/* Looks at all the particles of the parent box and distributes them amongst 
   the children. If the child does not exist, one is created. */
void
SubdivideBox (int my_id, box *b)
{
   particle *p;
   box *child;
   int i;

   for (i = 0; i < b->num_particles; i++) {
      p = b->particles[i];
      if (p->pos.y > b->y_center) {
	 if (p->pos.x > b->x_center) {
	    child = b->children[0];
	    if (child == NULL)
	       child = CreateChild(my_id, b, 0);
	 }
	 else {
	    child = b->children[1];
	    if (child == NULL)
	       child = CreateChild(my_id, b, 1);
	 }
      }
      else {
	 if (p->pos.x > b->x_center) {
	    child = b->children[3];
	    if (child == NULL)
	       child = CreateChild(my_id, b, 3);
	 }
	 else {
	    child = b->children[2];
	    if (child == NULL)
	       child = CreateChild(my_id, b, 2);
	 }
      }
      child->particles[child->num_particles++] = p;
      b->particles[i] = NULL;
   }
   b->num_particles = 0;
   b->type = PARENT;
}


/* Each processor keeps track of the boxes that it has inserted into the tree.
   This list is called a partition because when construction is over, every box
   in the global tree will also reside in one and only one of the processors'
   partitions. This is needed for list construction (which you don't have to 
   worry about) and cost zone computation (which you will). */
void
MergeLocalGrid (int my_id)
{
   MLGHelper(my_id, Local[my_id].Local_Grid, Grid, NULL);
}


void
MLGHelper (int my_id, box *local_box, box *global_box, box *global_parent)
{
   int success;
   int i;

   success = FALSE;
   while (success == FALSE) {
      if (local_box->type == PARENT) {
	 if (global_box == NULL) {
	    success = InsertBoxInGrid(my_id, local_box, global_parent);
	 }
	 else {
	    if (global_box->type == PARENT) {
	       success = TRUE;
	       for (i = 0; i < NUM_OFFSPRING; i++) {
		  if (local_box->children[i] != NULL)
		     MLGHelper(my_id, local_box->children[i],
			       global_box->children[i], global_box);
	       }
	    }
	    else {
	       success = RemoveBoxFromGrid(my_id, global_box, global_parent);
	       if (success == TRUE) {
		  InsertParticlesInTree(my_id, global_box->particles, 
					global_box->num_particles, local_box);
		  success = InsertBoxInGrid(my_id, local_box, global_parent);
	       }
	    }
	 }
      }
      else {
	 if (global_box == NULL) {
	    success = InsertBoxInGrid(my_id, local_box, global_parent);
	 }
	 else {
	    if (global_box->type == PARENT) {
	       MergeLocalParticles(my_id, local_box->particles,
		                   local_box->num_particles, global_box);
	       success = TRUE;
	    }
	    else {
	       success = RemoveBoxFromGrid(my_id, global_box, global_parent);
	       if (success == TRUE) {
		  InsertParticlesInLeaf(my_id, global_box->particles, 
					global_box->num_particles, local_box);
		  success = InsertBoxInGrid(my_id, local_box, global_parent);
	       }
	    }
	 }
      }
      if (success == FALSE) {
	 if (global_parent == NULL)
	    global_box = Grid;
	 else
	    global_box = global_parent->children[local_box->child_num];
      }
   }
}


void
MergeLocalParticles (int my_id, particle **p_array, int num_of_particles, box *pb)
{
   particle *(p_dist)[NUM_OFFSPRING][MAX_PARTICLES_PER_BOX];
   int num_p_dist[NUM_OFFSPRING];
   box *child;
   box *new_box;
   int success;
   int i, j, k;

   SplitParticles(my_id, p_array, num_of_particles, 
		  (particle **) p_dist, num_p_dist, pb);
   for (i= 0; i < NUM_OFFSPRING; i++) {
      if (num_p_dist[i] > 0) {
	 child = pb->children[i];
	 if (child == NULL) {
	    child = CreateLeaf(my_id, pb, i, p_dist[i], num_p_dist[i]);
	    success = InsertBoxInGrid(my_id, child, pb);
	 }
	 else {
	    if (child->type == PARENT) {
	       MergeLocalParticles(my_id, p_dist[i], num_p_dist[i], child);
	       success = TRUE;
	    }
	    else {
	       success = RemoveBoxFromGrid(my_id, child, pb);
	       if (success == TRUE) {
		  InsertParticlesInLeaf(my_id, p_dist[i], num_p_dist[i], child);
		  success = InsertBoxInGrid(my_id, child, pb);
	       }
	       else
		  child = CreateLeaf(my_id, pb, i, p_dist[i], num_p_dist[i]);
	    }
	 }
	 if (success == FALSE) {
	    MLGHelper(my_id, child, pb->children[child->child_num], pb);
	 }
      }
   }
}


void
SplitParticles (int my_id, particle **p_array, int length, particle **p_dist,
		int num_p_dist[NUM_OFFSPRING], box *pb)
{
   particle *p;
   int i;
  
   for (i = 0; i < NUM_OFFSPRING; i++)
      num_p_dist[i] = 0;
   for (i = 0; i < length; i++) {
      p = p_array[i];
      if (p->pos.y > pb->y_center) {
	 if (p->pos.x > pb->x_center)
	    *(p_dist + num_p_dist[0]++) = p;
	 else 
	    *(p_dist + MAX_PARTICLES_PER_BOX + num_p_dist[1]++) = p;
      }		
      else {
	 if (p->pos.x > pb->x_center)
	    *(p_dist + (3 * MAX_PARTICLES_PER_BOX) + num_p_dist[3]++) = p;
	 else 
	    *(p_dist + (2 * MAX_PARTICLES_PER_BOX) + num_p_dist[2]++) = p;
      }	  
   }
}


box *
CreateLeaf (int my_id, box *pb, int new_child_num, particle **p_array, int length)
{
   real child_length, child_offset;
   box *ret_box;
   int i;
    
   child_length = pb->length / (real) NUM_DIMENSIONS;
   child_offset = pb->length / (real) NUM_OFFSPRING;  
   if (new_child_num == 0) {
      ret_box = InitBox(my_id, (pb->x_center + child_offset), 
			(pb->y_center + child_offset), child_length, 
			pb);
   }
   if (new_child_num == 1) {
      ret_box = InitBox(my_id, (pb->x_center - child_offset), 
			(pb->y_center + child_offset), child_length, 
			pb);
   }
   if (new_child_num == 2) {
      ret_box = InitBox(my_id, (pb->x_center - child_offset), 
			(pb->y_center - child_offset), child_length, 
			pb);
   }
   if (new_child_num == 3) {
      ret_box = InitBox(my_id, (pb->x_center + child_offset), 
			(pb->y_center - child_offset), child_length, 
			pb);
   }
   ret_box->child_num = new_child_num;
   for (i = 0; i < length; i++) {
      ret_box->particles[i] = p_array[i];
   }
   ret_box->num_particles = length;
   return ret_box;
}


void
InsertParticlesInLeaf (int my_id, particle **p_array, int length, box *b)
{
   int i, j;
   int offset;
  
   if ((length + b->num_particles) > MAX_PARTICLES_PER_BOX) {
      for (i = b->num_particles, j = length - 1; i < MAX_PARTICLES_PER_BOX;
	   i++, j--)
	 b->particles[i] = p_array[j];
      b->num_particles = MAX_PARTICLES_PER_BOX;
      length = j + 1;
      InsertParticlesInTree(my_id, p_array, length, b);
   }
   else {
      offset = b->num_particles;
      for (i = 0; i < length; i++)
	 b->particles[i + offset] = p_array[i];
      b->num_particles += length;
   }
}


int
InsertBoxInGrid (int my_id, box *b, box *pb)
{
   int success;
  
   if (pb == NULL) {
      {pthread_mutex_lock(&(G_Memory->single_lock));};
      if (Grid == NULL) {
	 Grid = b;
	 success = TRUE;
      }
      else
	 success = FALSE;
      {pthread_mutex_unlock(&(G_Memory->single_lock));};
   }
   else {
      {pthread_mutex_lock(&((G_Memory->lock_array)[(pb->particle_lock_index)]));};
      if (pb->children[b->child_num] == NULL) {
	 pb->children[b->child_num] = b;
	 pb->num_children += 1;
	 b->parent = pb;
	 success = TRUE;
      }
      else
	 success = FALSE;
      {pthread_mutex_unlock(&((G_Memory->lock_array)[(pb->particle_lock_index)]));};
   }
   if (success == TRUE)
      InsertSubtreeInPartition(my_id, b);
   return success;
}


int
RemoveBoxFromGrid (int my_id, box *b, box *pb)
{
   int success;
  
   if (pb == NULL) {
      {pthread_mutex_lock(&(G_Memory->single_lock));};
      if (Grid == b) {
	 Grid = NULL;
	 success = TRUE;
      }
      else
	 success = FALSE;
      {pthread_mutex_unlock(&(G_Memory->single_lock));};
   }
   else {
      {pthread_mutex_lock(&((G_Memory->lock_array)[(pb->particle_lock_index)]));};
      if (pb->children[b->child_num] == b) {
	 pb->children[b->child_num] = NULL;
	 b->parent = NULL;
	 pb->num_children -= 1;
	 success = TRUE;
      }
      else
	 success = FALSE;
      {pthread_mutex_unlock(&((G_Memory->lock_array)[(pb->particle_lock_index)]));};
   }
   return success;
}


void
InsertSubtreeInPartition (int my_id, box *b)
{
   int i;
   box *child;
  
   if (b->proc == my_id) {
      InsertBoxInPartition(my_id, b);
   }
   if (b->type == PARENT) {
      for (i = 0; i < NUM_OFFSPRING; i++) {
	 child = b->children[i];
	 if (child == NULL)
	    child = b->shadow[i];
	 if (child != NULL)
	    InsertSubtreeInPartition(my_id, child);
      }
   }
}


void
CleanupGrid (int my_id)
{
   box *b_scan, *tb;
   int i;

   b_scan = Local[my_id].Childless_Partition;
   while (b_scan != NULL) {
      if (((b_scan->parent != NULL) || (b_scan == Grid))
	  && (b_scan->type == CHILDLESS))
	 b_scan = b_scan->next;
      else {
	 tb = b_scan;
	 b_scan = b_scan->next;
	 if (tb->type == PARENT) {
	    tb->type = CHILDLESS;
	    RemoveBoxFromPartition(my_id, tb);
	    tb->type = PARENT;
	    if ((tb->parent != NULL) || (tb == Grid)) {
	       InsertBoxInPartition(my_id, tb);
	    }
	 }
	 else
	    RemoveBoxFromPartition(my_id, tb);
      }
   }
}


void
ConstructGridLists (int my_id, box *b)
{
   SetSiblings(my_id, b);
   SetColleagues(my_id, b);
}


void
SetSiblings (int my_id, box *b)
{
   box *pb, *sb;
   int i;
  
   b->num_siblings = 0;
   pb = b->parent;
   if (pb != NULL) {
      for (i = 0; i < NUM_OFFSPRING; i++) {
	 sb = pb->children[i];
	 if ((sb != NULL) && (sb != b))
	    b->siblings[b->num_siblings++] = sb;
      }
   }
}


void
SetColleagues (int my_id, box *b)
{
   box *pb, *cb, *cousin;
   int i, j;

   b->num_colleagues = 0;
   pb = b->parent;
   if (pb != NULL) {
      for (i = 0; i < b->num_siblings; i++)
	 b->colleagues[b->num_colleagues++] = b->siblings[i];
      while (b->construct_synch == 0) {
	 /* wait */;
      }
      b->construct_synch = 0;
      for (i = 0; i < pb->num_colleagues; i++) {
	 cb = pb->colleagues[i];
	 for (j = 0; j < NUM_OFFSPRING; j++) {
	    cousin = cb->children[j];
	    if (cousin != NULL) {
	       if (AdjacentBoxes(my_id, b, cousin) == TRUE)
		  b->colleagues[b->num_colleagues++] = cousin;
	    }
	 }
      }
   }
   if (b->type == PARENT) {
      for (i = 0; i < NUM_OFFSPRING; i++) {
	 if (b->children[i] != NULL) {
	    b->children[i]->construct_synch = 1;
	 }
      }
   }
}


/*
 *  ConstructInteractionLists (int my_id, box *b)
 *
 *  Args : a box, b.
 *
 *  Returns : nothing.
 *
 *  Side Effects : Creates b's colleagues, u_list,
 *     v_list, and w_list.
 *
 */
void
ConstructInteractionLists (int my_id, box *b)
{

   SetVList(my_id, b);
   if (b->type == CHILDLESS) {
      SetUList(my_id, b);
      SetWList(my_id, b);
   }

}


void
SetVList (int my_id, box *b)
{
   box *pb, *cb, *cousin;
   int i, j;
  
   b->num_v_list = 0;
   pb = b->parent;
   if (pb != NULL) {
      for (i = 0; i < pb->num_colleagues; i++) {
	 cb = pb->colleagues[i];
	 for (j = 0; j < NUM_OFFSPRING; j++) {
	    cousin = cb->children[j];
	    if (cousin != NULL) {
	       if (WellSeparatedBoxes(my_id, b, cousin) == TRUE)
		  b->v_list[b->num_v_list++] = cousin;
	    }
	 }
      }
   }
}  


/*
 *  SetUList (int my_id, box *b)
 *
 *  Args : a box, b.
 *
 *  Returns : nothing.
 *
 *  Side Effects : Sets b's adjacency interaction list.
 *
 *  Comments : The helper function is used to enable recursion.
 *
 */
void
SetUList (int my_id, box *b)
{
   b->num_u_list = 0;
   SetUListHelper(my_id, b, Grid);

}


/*
 *  SetUListHelper (int my_id, box *b, box *pb)
 *
 *  Args : a box, b, and a parent box, pb.
 *
 *  Returns : nothing.
 *
 *  Side Effects : Sets b's adjacency interaction list.
 *
 *  Comments :
 *
 */
void
SetUListHelper (int my_id, box *b, box *pb)
{
   box *child;
   int i;

   for (i = 0; i < NUM_OFFSPRING; i++) {
      child = pb->children[i];
      if (child != NULL) {
	 if (AdjacentBoxes(my_id, b, child) == TRUE) {
	    if (child->type == CHILDLESS)
	       b->u_list[b->num_u_list++] = child;
	    else
	       SetUListHelper(my_id, b, child);
	 }
	 else {
	    if (AncestorBox(my_id, b, child) == TRUE)
	       SetUListHelper(my_id, b, child);
	 }
      }
   }

}


/*
 *  AncestorBox (int my_id, box *b, box *ancestor_box)
 *
 *  Args : a box, b, and its possible ancestor box, ancestor_box.
 *
 *  Returns : TRUE, if ancestor_box is an ancestor of b, FALSE if not.
 *
 *  Side Effects : none.
 *
 *  Comments : A box is NOT the ancestor of himself. So, AncestorBox(my_id, b,b)
 *    always returns FALSE. So, ancestor_box is indeed the ancestor of b
 *    if their sizes are not equal and if b's center lies within the boundaries
 *    of ancestor_box.
 *
 */
int
AncestorBox (int my_id, box *b, box *ancestor_box)
{
   real x_center_distance;
   real y_center_distance;
   int ret_val = TRUE;

   if (b->length != ancestor_box->length) {

      x_center_distance = fabs((double)(b->x_center - ancestor_box->x_center));
      y_center_distance = fabs((double)(b->y_center - ancestor_box->y_center));
      if ((x_center_distance > (ancestor_box->length / 2.0)) ||
	  (y_center_distance > (ancestor_box->length / 2.0)))
	 ret_val = FALSE;
   }
   else
      ret_val = FALSE;

   return ret_val;

}


/*
 *  SetWList (int my_id, box *b)
 *
 *  Args : a box, b.
 *
 *  Returns : nothing.
 *
 *  Side Effects : Sets b's w_list.
 *
 *  Comments : This list contains all descedants of b's colleagues whose
 *     parents are adjacent to b, but who are not adjacent to b themeselves.
 *     This list should only be computed for childless boxes.
 *
 */
void
SetWList (int my_id, box *b)
{
   box *co_search;
   int i;

   b->num_w_list = 0;
   for (i = 0; i < b->num_colleagues; i++) {
      co_search = b->colleagues[i];
      if (co_search->type == PARENT)
	 InsertNonAdjChildren(my_id, b, co_search);
   }

}


/*
 *  InsertNonAdjChildren (int my_id, box *b, box *pb)
 *
 *  Args : a parent box, pb, and the box with the weak iteraction list, b.
 *
 *  Returns : nothing.
 *
 *  Side Effects : Inserts all non-adjacent children of adjacent box pb into
 *     b's weak interaction list.
 *
 *  Comments : This function recursively calls itself on every child that is
 *     also adjacent to b and a parent.
 *
 */
void
InsertNonAdjChildren (int my_id, box *b, box *pb)
{
   int i;
   box *child;

   for (i = 0; i < pb->num_children; i++) {
      child = pb->children[i];
      if (child != NULL) {
	 if (AdjacentBoxes(my_id, b, child) == TRUE) {
	    if (child->type == PARENT)
	       InsertNonAdjChildren(my_id, b, child);
	 }
	 else
	    b->w_list[b->num_w_list++] = child;
      }
   }
  
}



/* Generated from ../Source/construct_grid.C */
