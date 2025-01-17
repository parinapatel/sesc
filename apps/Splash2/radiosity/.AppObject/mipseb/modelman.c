

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

/************************************************************************
 *
 *	Model definition & management.
 *
 *       This module has the following functions:
 *       (1) Create modeling tasks.
 *       (2) Process modeling tasks.
 *
 *************************************************************************/

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


/***************************************************************************
 *
 *    Test Model definitions.
 *
 ****************************************************************************/

#define  WHITE  { 0.9, 0.9, 0.9 }
#define  GREEN  { 0.1, 0.9, 0.1 }
#define  PASSIVE { 0.0, 0.0, 0.0 }
#define  WHITE_LIGHT { 200.0, 200.0, 200.0 }

Model_Triangle  tri1 = { WHITE,
                             PASSIVE,
                             { 0.0, 0.0, 0.0 },
                             { 100.0, 0.0, 0.0 },
                             { 0.0, 100.0, 0.0 } } ;

Model_Triangle  tri2 = { WHITE,
                             PASSIVE,
                             { 0.0, 0.0, 100.0 },
                             { 100.0, 100.0, 0.0 },
                             { 0.0, 0.0, -100.0 } } ;

Model_Rectangle  Floor = { WHITE,
                               PASSIVE,
                               {  600.0, -600.0,  700.0 },
                               {  600.0, -600.0, -700.0 },
                               { -800.0, -600.0,  700.0 } } ;

Model_Rectangle  Ceiling = { WHITE,
                                 PASSIVE,
                                 {  600.0,  600.0,  700.0 },
                                 { -800.0,  600.0,  700.0 },
                                 {  600.0,  600.0, -700.0 } } ;

Model_Rectangle  wall1 = { WHITE,
                               PASSIVE,
                               { -800.0, -600.0,  700.0 },
                               { -800.0, -600.0, -700.0 },
                               { -800.0,  600.0,  700.0 } } ;

Model_Rectangle  wall2 = { WHITE,
                               PASSIVE,
                               {  600.0, -600.0, -700.0 },
                               {  600.0,  600.0, -700.0 },
                               { -800.0, -600.0, -700.0 } } ;

Model_Rectangle  deskTop =  { GREEN,
                                  PASSIVE,
                                  { -795.0, -320.0,  300.0 },
                                  { -400.0, -320.0,  300.0 },
                                  { -795.0, -320.0, -300.0 } } ;

Model_Rectangle  deskBtm =  { GREEN,
                                  PASSIVE,
                                  { -795.0, -340.0,  300.0 },
                                  { -795.0, -340.0, -300.0 },
                                  { -400.0, -340.0,  300.0 } } ;

Model_Rectangle  light1 = { WHITE,
                                WHITE_LIGHT,
                                { -795.0, -50.0,  50.0 },
                                { -795.0, -50.0, -50.0 },
                                { -795.0,  50.0,  50.0 } } ;


/***************************************************************************
 ****************************************************************************
 *
 *    Methods of Model object
 *
 ****************************************************************************
 ****************************************************************************/
/***************************************************************************
 *
 *    init_modeling_tasks()
 *
 *    Initialize modeling task.
 *
 ****************************************************************************/

int model_selector = MODEL_TEST_DATA ;
static void init_test_model_tasks(), init_room_model_tasks() ;
  static void init_largeroom_model_tasks() ;
  
  void init_modeling_tasks(process_id)
  unsigned process_id;
{
    
    extern ModelDataBase room_model[] ;
    extern ModelDataBase largeroom_model[] ;
    
    if( ! check_task_counter(process_id) )
        return ;
    
    switch( model_selector )
        {
        case MODEL_TEST_DATA:
        default:
            init_test_model_tasks(process_id) ;
            break ;
        case MODEL_ROOM_DATA:
            init_room_model_tasks( room_model, process_id ) ;
            break ;
        case MODEL_LARGEROOM_DATA:
            init_room_model_tasks( largeroom_model, process_id ) ;
            break ;
        }
}



static void init_test_model_tasks(process_id)
  unsigned process_id;
{
    create_modeling_task( (Model*)&Floor, MODEL_RECTANGLE, process_id ) ;
    create_modeling_task( (Model*)&Ceiling, MODEL_RECTANGLE, process_id ) ;
    create_modeling_task( (Model*)&wall1, MODEL_RECTANGLE, process_id ) ;
    create_modeling_task( (Model*)&wall2, MODEL_RECTANGLE, process_id ) ;
    create_modeling_task( (Model*)&deskTop, MODEL_RECTANGLE, process_id ) ;
    create_modeling_task( (Model*)&deskBtm, MODEL_RECTANGLE, process_id ) ;
    create_modeling_task( (Model*)&light1, MODEL_RECTANGLE, process_id ) ;
}



static void init_room_model_tasks( model, process_id )
  
  ModelDataBase *model ;
  unsigned process_id;
{
    ModelDataBase *pm ;
    
    for( pm = model ; pm->type != MODEL_NULL ; pm++ )
        create_modeling_task( &pm->model, pm->type, process_id ) ;
}


#define MODEL_START (140)

static void init_room_model_tasks2( model, process_id )
  
  ModelDataBase *model ;
  unsigned process_id;
{
    ModelDataBase *pm ;
    int cnt = 0 ;
    pm = &model[ cnt ] ;
    
    for( ; pm->type != MODEL_NULL ; pm++ )
        {
            if( pm->model.emittance.g == 0 )
                pm->model.emittance.g = 0.5 + (cnt & 3)*0.1  ;
            else
                pm->model.emittance.g = 1.0 ;
            if( cnt >= MODEL_START )
                create_modeling_task( &pm->model, pm->type, process_id ) ;
            if( pm->type == MODEL_TRIANGLE )
                {
                    pm++ ;
                    if( pm->model.emittance.g == 0 )
                        pm->model.emittance.g = 0.5 + (cnt & 3)*0.1  ;
                    else
                        pm->model.emittance.g = 1.0 ;
                    if( cnt >= MODEL_START )
                        create_modeling_task( &pm->model, pm->type, process_id ) ;
                }
            
            if( cnt >= MODEL_START )
                {
                    printf( "Polygon %d\n", cnt ) ;
                    g_clear() ;
                    display_scene( 1, 1, 0, 0, process_id ) ;
                    while( getchar() != '\n' ) ;
                }
            cnt++ ;
        }
}



/***************************************************************************
 *
 *    process_model()
 *
 ****************************************************************************/

static void process_triangle(), process_rectangle() ;
  
  
  void process_model( model, type, process_id )
  
  Model *model ;
  int    type ;
  unsigned process_id;
{
    switch( type )
        {
        case MODEL_TRIANGLE:
            process_triangle( (Model_Triangle *)model, process_id ) ;
            break ;
        case MODEL_RECTANGLE:
            process_rectangle( (Model_Rectangle *)model, process_id ) ;
            break ;
        default:
            fprintf( stderr, "Panic:process_model:Illegal type %d\n", type ) ;
        }
}



/***************************************************************************
 *
 *    process_triangle()
 *
 ****************************************************************************/

static void process_triangle( model, process_id )
  
  Model_Triangle *model ;
  unsigned process_id;
{
    Patch *p ;
    float length ;
    
    /* Create a patch */
    p = get_patch(process_id) ;
    
    /* (1) Set the Vertecies */
    p->p1  = model->p1 ;
    p->p2  = model->p2 ;
    p->p3  = model->p3 ;
    
    /* (2) Create the Edges */
    p->ev1 = create_elemvertex( &p->p1, process_id ) ;
    p->ev2 = create_elemvertex( &p->p2, process_id ) ;
    p->ev3 = create_elemvertex( &p->p3, process_id ) ;
    p->e12 = create_edge( p->ev1, p->ev2, process_id ) ;
    p->e23 = create_edge( p->ev2, p->ev3, process_id ) ;
    p->e31 = create_edge( p->ev3, p->ev1, process_id ) ;
    
    /* (3) Other patch properties */
    length = comp_plane_equ( &p->plane_equ,
                            &model->p1, &model->p2, &model->p3, process_id ) ;
    p->area      = length * (float)0.5 ;
    p->color     = model->color ;
    p->emittance = model->emittance ;
    
    
    /* Create a BSP insertion task */
    create_bsp_task( p, global->bsp_root, process_id ) ;
}


/***************************************************************************
 *
 *    process_rectangle()
 *
 ****************************************************************************/

static void process_rectangle( model, process_id )
  
  Model_Rectangle *model ;
  unsigned process_id;
{
    Patch *p, *q ;
    float length ;
    
    /* Create a patch (P1-P2-P3) */
    p = get_patch(process_id) ;
    
    /* (1) Set the Vertecies */
    p->p1  = model->p1 ;
    p->p2  = model->p2 ;
    p->p3  = model->p3 ;
    
    /* (2) Create the Edges */
    p->ev1 = create_elemvertex( &p->p1, process_id ) ;
    p->ev2 = create_elemvertex( &p->p2, process_id ) ;
    p->ev3 = create_elemvertex( &p->p3, process_id ) ;
    p->e12 = create_edge( p->ev1, p->ev2, process_id ) ;
    p->e23 = create_edge( p->ev2, p->ev3, process_id ) ;
    p->e31 = create_edge( p->ev3, p->ev1, process_id ) ;
    
    /* (3) Other patch properties */
    length = comp_plane_equ( &p->plane_equ,
                            &model->p1, &model->p2, &model->p3, process_id ) ;
    p->area     = length * (float)0.5 ;
    p->color    = model->color ;
    p->emittance= model->emittance ;
    
    
    /* Create a patch (P(2+3-1)-P3-P2) */
    q = get_patch(process_id) ;
    
    /* (1) Set the Vertices */
    q->p1.x     = model->p2.x + model->p3.x - model->p1.x ;
    q->p1.y     = model->p2.y + model->p3.y - model->p1.y ;
    q->p1.z     = model->p2.z + model->p3.z - model->p1.z ;
    q->p2       = model->p3 ;
    q->p3       = model->p2 ;
    
    /* (2) Create the Edges */
    q->ev1 = create_elemvertex( &q->p1, process_id ) ;
    q->ev2 = p->ev3 ;
    q->ev3 = p->ev2 ;
    q->e12 = create_edge( q->ev1, q->ev2, process_id ) ;
    q->e23 = p->e23 ;
    q->e31 = create_edge( q->ev3, q->ev1, process_id ) ;
    
    /* (3) Other patch properties */
    q->plane_equ= p->plane_equ ;
    q->area     = p->area ;
    q->color    = p->color ;
    q->emittance= p->emittance ;
    
    /* Create BSP insertion tasks */
    create_bsp_task( p, global->bsp_root, process_id ) ;
    create_bsp_task( q, global->bsp_root, process_id ) ;
}



/* Generated from ../Changed/modelman.C */
