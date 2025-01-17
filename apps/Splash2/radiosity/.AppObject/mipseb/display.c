

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

#include <stdint.h>
#include <stdio.h>



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
;

#include "radiosity.h"
  
  
/************************************************************************
 *
 *    radiosity_averaging
 *
 *************************************************************************/
  
  static void add_radiosity_to_vertex() ;

void radiosity_averaging( elem, mode, process_id )
  
  Element *elem ;
  int mode ;
  unsigned process_id;
{
    float inv_weight ;
    Vertex pc ;
    int reverse ;
    
    if( ! LEAF_ELEMENT(elem) )
        {
            create_radavg_task( elem->center, mode, process_id ) ;
            create_radavg_task( elem->top,    mode, process_id ) ;
            create_radavg_task( elem->right,  mode, process_id ) ;
            create_radavg_task( elem->left,   mode, process_id ) ;
            return ;
        }
    
    else if( mode == RAD_AVERAGING_MODE )
        {
            /* Compute center point */
            center_point( &elem->ev1->p, &elem->ev2->p, &elem->ev3->p, &pc ) ;
            
            reverse = EDGE_REVERSE( elem->e12, elem->ev1, elem->ev2 ) ;
            foreach_leaf_edge( elem->e12, reverse,
                              add_radiosity_to_vertex, (uintptr_t)elem, (uintptr_t)&pc, process_id ) ;
            reverse = EDGE_REVERSE( elem->e23, elem->ev2, elem->ev3 ) ;
            foreach_leaf_edge( elem->e23, reverse,
                              add_radiosity_to_vertex, (uintptr_t)elem, (uintptr_t)&pc, process_id ) ;
            reverse = EDGE_REVERSE( elem->e31, elem->ev3, elem->ev1 ) ;
            foreach_leaf_edge( elem->e31, reverse,
                              add_radiosity_to_vertex, (uintptr_t)elem, (uintptr_t)&pc, process_id ) ;
        }
    else
        {
            /* Normalize it */
            {pthread_mutex_lock(&(elem->ev1->ev_lock->lock));};
            if( elem->ev1->weight != 1.0 )
                {
                    inv_weight = (float)1.0 / elem->ev1->weight  ;
                    elem->ev1->col.r *= inv_weight ;
                    elem->ev1->col.g *= inv_weight ;
                    elem->ev1->col.b *= inv_weight ;
                    elem->ev1->weight = 1.0 ;
                }
            {pthread_mutex_unlock(&(elem->ev1->ev_lock->lock));};
            
            {pthread_mutex_lock(&(elem->ev2->ev_lock->lock));};
            if( elem->ev2->weight != 1.0 )
                {
                    inv_weight = (float)1.0 / elem->ev2->weight  ;
                    elem->ev2->col.r *= inv_weight ;
                    elem->ev2->col.g *= inv_weight ;
                    elem->ev2->col.b *= inv_weight ;
                    elem->ev2->weight = 1.0 ;
                }
            {pthread_mutex_unlock(&(elem->ev2->ev_lock->lock));};
            
            {pthread_mutex_lock(&(elem->ev3->ev_lock->lock));};
            if( elem->ev3->weight != 1.0 )
                {
                    inv_weight = (float)1.0 / elem->ev3->weight  ;
                    elem->ev3->col.r *= inv_weight ;
                    elem->ev3->col.g *= inv_weight ;
                    elem->ev3->col.b *= inv_weight ;
                    elem->ev3->weight = 1.0 ;
                }
            {pthread_mutex_unlock(&(elem->ev3->ev_lock->lock));};
        }
}

static void add_radiosity_to_vertex( edge, reverse, elem, p_c, process_id )
  
  Edge *edge ;
  int reverse ;			/* Direction */
  Element *elem ;
  Vertex *p_c ;
  unsigned process_id;
{
    ElemVertex *ev ;
    float weight ;
    
    if( reverse )
        ev = edge->pb ;
    else
        ev = edge->pa ;
    
    weight = (float)1.0 / distance( &ev->p, p_c ) ;
    weight = 1.0 ;
    weight = elem->area ;
    {pthread_mutex_lock(&(ev->ev_lock->lock));};
    ev->col.r += (elem->rad.r * weight) ;
    ev->col.g += (elem->rad.g * weight) ;
    ev->col.b += (elem->rad.b * weight) ;
    ev->weight += weight ;
    {pthread_mutex_unlock(&(ev->ev_lock->lock));};
}


/************************************************************************
 *
 *    setup_view()
 *
 *************************************************************************/

Vertex view_vec ;  /* Origin to the viewer */
static float view_rot_x, view_rot_y, view_dist, view_zoom ;


void setup_view( rot_x, rot_y, dist, zoom, process_id )
  
  float rot_x, rot_y, dist, zoom ;
  unsigned process_id;
{
    Vertex v1, v2 ;
    float cc, ss ;
    
    /* Save parameters */
    view_rot_x = rot_x ;
    view_rot_y = rot_y ;
    view_dist  = dist ;
    view_zoom  = zoom ;
    
    /* Compute view vector */
    v1.x = 0.0 ;
    v1.y = 0.0 ;
    v1.z = 1.0 ;
    
    /* Rotate view vector */
    cc = cos( -rot_x * (M_PI / 180.0) ) ;
    ss = sin( -rot_x * (M_PI / 180.0) ) ;
    v2.x = v1.x ;
    v2.y = cc * v1.y - ss * v1.z ;
    v2.z = ss * v1.y + cc * v1.z ;
    
    cc = cos( -rot_y * (M_PI / 180.0) ) ;
    ss = sin( -rot_y * (M_PI / 180.0) ) ;
    v1.z = cc * v2.z - ss * v2.x ;
    v1.x = ss * v2.z + cc * v2.x ;
    v1.y = v2.y ;
    
    /* Store view vector */
    view_vec = v1 ;
}


/************************************************************************
 *
 *    display_scene()
 *
 *************************************************************************/

void display_scene( fill_sw, patch_sw, mesh_sw, interaction_sw, process_id )
  
  int fill_sw, patch_sw, mesh_sw, interaction_sw ;
  unsigned process_id;
{
    /* Clear the screen */
    g_clear() ;
    
    /* Set matrix */ 
    g_setup_view( view_rot_x, view_rot_y, view_dist, view_zoom ) ;
    
    if( fill_sw == 2 )
        {
            /* Fill surfaces */
            display_elements_in_bsp_tree( DISPLAY_SHADED, process_id ) ;
        }
    if( fill_sw == 1 )
        {
            /* Fill surfaces */
            display_elements_in_bsp_tree( DISPLAY_FILLED, process_id ) ;
        }
    if( mesh_sw )
        {
            /* Draw mesh */
            g_color( G_BLUE ) ;
            display_elements_in_bsp_tree( DISPLAY_EDGEONLY, process_id ) ;
        }
    if( patch_sw )
        {
            g_color( G_RED ) ;
            display_patches_in_bsp_tree( DISPLAY_EDGEONLY, process_id ) ;
        }
    if( interaction_sw )
        {
            g_color( G_GREEN ) ;
            display_interactions_in_bsp_tree(process_id) ;
        }
    
    /* Flush */
    g_flush() ;
}

/************************************************************************
 *
 *    display_patch()
 *
 *************************************************************************/

void display_patch( patch, mode, process_id )
  
  Patch *patch ;
  int   mode ;
  unsigned process_id;
{
    Vertex p_buf[4] ;
    Rgb   c_buf[4] ;
    
    if( mode == DISPLAY_SHADED )
        {
            if( inner_product( &patch->plane_equ.n, &view_vec ) < F_ZERO )
                return ;
            
            p_buf[0] = patch->p1 ;
            p_buf[1] = patch->p2 ;
            p_buf[2] = patch->p3 ;
            c_buf[0] = patch->color ;
            c_buf[1] = patch->color ;
            c_buf[2] = patch->color ;
            
            g_spolygon( 3, p_buf, c_buf ) ;
        }
    else if( mode == DISPLAY_FILLED )
        {
            if( inner_product( &patch->plane_equ.n, &view_vec ) < F_ZERO )
                return ;
            
            p_buf[0] = patch->p1 ;
            p_buf[1] = patch->p2 ;
            p_buf[2] = patch->p3 ;
            
            g_polygon( 3, p_buf ) ;
        }
    else
        {
            g_line( &patch->p1, &patch->p2 ) ;
            g_line( &patch->p2, &patch->p3 ) ;
            g_line( &patch->p3, &patch->p1 ) ;
        }
}


/************************************************************************
 *
 *    display_patches_in_bsp_tree()
 *
 *************************************************************************/

void display_patches_in_bsp_tree( mode, process_id )
  
  int mode ;
  unsigned process_id;
{
    foreach_depth_sorted_patch( &view_vec, display_patch, (int)mode, process_id ) ;
}



/************************************************************************
 *
 *    display_element()
 *
 *************************************************************************/

static void _display_shaded_triangle() ;
  
  void display_element( element, mode, process_id )
  
  Element *element ;
  int   mode ;
  unsigned process_id;
{
    Vertex p_buf[4] ;
    
    if( inner_product( &element->patch->plane_equ.n, &view_vec ) < F_ZERO )
        return ;
    
    if( mode == DISPLAY_SHADED )
        {
            _display_shaded_triangle( element->ev1, element->ev2,
                                     element->ev3,
                                     element->e12, element->e23, element->e31, process_id ) ;
        }
    else if( mode == DISPLAY_FILLED )
        {
            g_rgb( element->rad ) ;
            p_buf[0] = element->ev1->p ;
            p_buf[1] = element->ev2->p ;
            p_buf[2] = element->ev3->p ;
            
            g_polygon( 3, p_buf ) ;
        }
    else
        {
            g_line( &element->ev1->p, &element->ev2->p ) ;
            g_line( &element->ev2->p, &element->ev3->p ) ;
            g_line( &element->ev3->p, &element->ev1->p ) ;
        }
}

static void _display_shaded_triangle( ev1, ev2, ev3, e12, e23, e31, process_id )
  
  ElemVertex *ev1, *ev2, *ev3 ;
  Edge *e12, *e23, *e31 ;
  unsigned process_id;
{
    Vertex p_buf[4] ;
    Rgb   c_buf[4] ;
    
    p_buf[0] = ev1->p ;
    p_buf[1] = ev2->p ;
    p_buf[2] = ev3->p ;
    c_buf[0] = ev1->col ;
    c_buf[1] = ev2->col ;
    c_buf[2] = ev3->col ;
    g_spolygon( 3, p_buf, c_buf ) ;
}


/************************************************************************
 *
 *    display_elements_in_patch()
 *
 *************************************************************************/

void display_elements_in_patch( patch, mode, process_id )
  
  Patch *patch ;
  int mode ;
  unsigned process_id;
{
    foreach_leaf_element_in_patch( patch, display_element, mode, process_id ) ;
    g_flush() ;
}


/************************************************************************
 *
 *    display_elements_in_bsp_tree()
 *
 *************************************************************************/

void display_elements_in_bsp_tree( mode, process_id )
  
  int mode ;
  unsigned process_id;
{
    foreach_depth_sorted_patch( &view_vec, display_elements_in_patch, mode, process_id );
}

/************************************************************************
 *
 *    display_interactions_in_element()
 *
 *************************************************************************/

static void _disp_interactions() ;
  
  void display_interactions_in_element( elem, mode, process_id )
  
  Element *elem ;
  int mode ;
  unsigned process_id;
{
    
    foreach_interaction_in_element( elem, _disp_interactions, mode, process_id ) ;
    g_flush() ;
}


static void _disp_interactions( elem, inter, mode, process_id )
  
  Element *elem ;
  Interaction *inter ;
  int mode ;
  unsigned process_id;
{
    Vertex pa, pb ;
    Element *edst ;
    
    
    /* Display interactions only with a particular patch */
    if(   (mode == DISPLAY_HALF_INTERACTIONS)
       && (inter->destination->patch->seq_no >= elem->patch->seq_no ) )
        return ;
    
    /* Compute mid point of the element */
    edst = inter->destination ;
    center_point( &elem->ev1->p, &elem->ev2->p, &elem->ev3->p, &pa ) ;
    center_point( &edst->ev1->p, &edst->ev2->p, &edst->ev3->p, &pb ) ;
    
    /* Draw a line */
    g_line( &pa, &pb ) ;
}



/************************************************************************
 *
 *    display_interactions_in_patch
 *
 *************************************************************************/

void display_interactions_in_patch( patch, mode, process_id )
  
  Patch *patch ;
  int mode ;  
  unsigned process_id;
{
    foreach_element_in_patch( patch, display_interactions_in_element, mode, process_id );
}

/************************************************************************
 *
 *    display_interactions_in_bsp_tree
 *
 *************************************************************************/

void display_interactions_in_bsp_tree(process_id)
  unsigned process_id;
{
    foreach_patch_in_bsp( display_interactions_in_patch,
                         DISPLAY_ALL_INTERACTIONS, process_id ) ;
}



/************************************************************************
 *************************************************************************
 *
 *  PostScript Version driver
 *
 *************************************************************************
 *************************************************************************/

/************************************************************************
 *
 *    ps_display_scene()
 *
 *************************************************************************/


void ps_display_scene( fill_sw, patch_sw, mesh_sw, interaction_sw, process_id )
  
  int fill_sw, patch_sw, mesh_sw, interaction_sw ;
  unsigned process_id;
{
    if( fill_sw )
        {
            /* Fill surfaces */
            ps_display_elements_in_bsp_tree( DISPLAY_SHADED, process_id ) ;
        }
    if( mesh_sw )
        {
            /* Draw mesh */
            ps_linewidth( 0.5 ) ;
            ps_display_elements_in_bsp_tree( DISPLAY_EDGEONLY, process_id ) ;
        }
    if( patch_sw )
        {
            /* Draw patches */
            ps_linewidth( 1.2 ) ;
            ps_display_patches_in_bsp_tree( DISPLAY_EDGEONLY, process_id ) ;
        }
    if( interaction_sw )
        {
            /* Draw interactions */
            ps_linewidth( 0.2 ) ;
            ps_display_interactions_in_bsp_tree(process_id) ;
        }
    
}

/************************************************************************
 *
 *    ps_display_patch()
 *
 *************************************************************************/

void ps_display_patch( patch, mode, process_id )
  
  Patch *patch ;
  int   mode ;
  unsigned process_id;
{
    Vertex p_buf[4] ;
    Rgb   c_buf[4] ;
    
    if( mode == DISPLAY_SHADED )
        {
            if( inner_product( &patch->plane_equ.n, &view_vec ) < F_ZERO )
                return ;
            p_buf[0] = patch->p1 ;
            p_buf[1] = patch->p2 ;
            p_buf[2] = patch->p3 ;
            c_buf[0] = patch->color ;
            c_buf[1] = patch->color ;
            c_buf[2] = patch->color ;
            
            ps_spolygon( 3, p_buf, c_buf ) ;
        }
    else if( mode == DISPLAY_FILLED )
        {
            if( inner_product( &patch->plane_equ.n, &view_vec ) < F_ZERO )
                return ;
            p_buf[0] = patch->p1 ;
            p_buf[1] = patch->p2 ;
            p_buf[2] = patch->p3 ;
            
            ps_polygon( 3, p_buf ) ;
        }
    else
        {
            p_buf[0] = patch->p1 ;
            p_buf[1] = patch->p2 ;
            p_buf[2] = patch->p3 ;
            
            ps_polygonedge( 3, p_buf ) ;
        }
}


/************************************************************************
 *
 *    ps_display_patches_in_bsp_tree()
 *
 *************************************************************************/

void ps_display_patches_in_bsp_tree( mode, process_id )
  
  int mode ;
  unsigned process_id;
{
    foreach_depth_sorted_patch( &view_vec, ps_display_patch, (int)mode, process_id ) ;
}



/************************************************************************
 *
 *    ps_display_element()
 *
 *************************************************************************/

void ps_display_element( element, mode, process_id )
  
  Element *element ;
  int   mode ;
  unsigned process_id;
{
    Vertex p_buf[4] ;
    Rgb   c_buf[4] ;
    
    if( mode == DISPLAY_SHADED )
        {
            if( inner_product( &element->patch->plane_equ.n, &view_vec )
               < F_ZERO )
                return ;
            p_buf[0] = element->ev1->p ;
            p_buf[1] = element->ev2->p ;
            p_buf[2] = element->ev3->p ;
            c_buf[0] = element->rad ;
            c_buf[1] = element->rad ;
            c_buf[2] = element->rad ;
            
            ps_spolygon( 3, p_buf, c_buf ) ;
        }
    else if( mode == DISPLAY_FILLED )
        {
            if( inner_product( &element->patch->plane_equ.n, &view_vec )
               < F_ZERO )
                return ;
            p_buf[0] = element->ev1->p ;
            p_buf[1] = element->ev2->p ;
            p_buf[2] = element->ev3->p ;
            
            ps_polygon( 3, p_buf ) ;
        }
    else
        {
            p_buf[0] = element->ev1->p ;
            p_buf[1] = element->ev2->p ;
            p_buf[2] = element->ev3->p ;
            
            ps_polygonedge( 3, p_buf ) ;
        }
}


/************************************************************************
 *
 *    ps_display_elements_in_patch()
 *
 *************************************************************************/

void ps_display_elements_in_patch( patch, mode, process_id )
  
  Patch *patch ;
  int mode ;
  unsigned process_id;
{
    foreach_leaf_element_in_patch( patch, ps_display_element, mode, process_id ) ;
}


/************************************************************************
 *
 *    ps_display_elements_in_bsp_tree()
 *
 *************************************************************************/

void ps_display_elements_in_bsp_tree( mode, process_id )
  
  int mode ;
  unsigned process_id;
{
    foreach_depth_sorted_patch( &view_vec, 
                               ps_display_elements_in_patch, mode, process_id ) ;
}

/************************************************************************
 *
 *    ps_display_interactions_in_element()
 *
 *************************************************************************/

static void _ps_disp_interactions() ;
  
  
  void ps_display_interactions_in_element( elem, mode, process_id )
  
  Element *elem ;
  int mode ;
  unsigned process_id;
{
    foreach_interaction_in_element( elem, _ps_disp_interactions, mode, process_id ) ;
}


static void _ps_disp_interactions( elem, inter, mode, process_id )
  
  Element *elem ;
  Interaction *inter ;
  int mode ;
  unsigned process_id;
{
    Vertex pa, pb ;
    Element *edst ;
    
    /* Display interactions only with a particular patch */
    if(   (mode == DISPLAY_HALF_INTERACTIONS)
       && (inter->destination->patch->seq_no >= elem->patch->seq_no ) )
        return ;
    
    /* Compute mid point of the element */
    edst = inter->destination ;
    center_point( &elem->ev1->p, &elem->ev2->p, &elem->ev3->p, &pa ) ;
    center_point( &edst->ev1->p, &edst->ev2->p, &edst->ev3->p, &pb ) ;
    
    /* Draw a line */
    ps_line( &pa, &pb ) ;
}



/************************************************************************
 *
 *    ps_display_interactions_in_patch
 *
 *************************************************************************/

void ps_display_interactions_in_patch( patch, mode, process_id )
  
  Patch *patch ;
  int mode ;  
  unsigned process_id;
{
    foreach_element_in_patch( patch,
                             ps_display_interactions_in_element, mode, process_id );
}

/************************************************************************
 *
 *    ps_display_interactions_in_bsp_tree
 *
 *************************************************************************/

void ps_display_interactions_in_bsp_tree(process_id)
  unsigned process_id;
{
    foreach_patch_in_bsp( ps_display_interactions_in_patch,
                         DISPLAY_ALL_INTERACTIONS, process_id ) ;
}


/* Generated from ../Changed/display.C */
