

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

/**************************************************************
 *
 *       Parallel Hierarchical Radiosity
 *
 *	    Main program
 *
 ***************************************************************/

#include <stdio.h>
#if defined(SGI_GL)
#include <gl.h>
#if defined(GL_NASA)
#include <panel.h>
#endif
#endif

/* ANL macro initialization */



#include <pthread.h>
#include <stdlib.h>
#include <semaphore.h>
#include <assert.h>
#if !(defined PAGE_SIZE)
#define PAGE_SIZE 4096
#endif
#define __MAX_THREADS__ 256

pthread_t __tid__[__MAX_THREADS__];
unsigned __threads__=0;
pthread_mutex_t __intern__;
void *our_malloc(size_t size, char * file, unsigned line) { return malloc(size); }
;

#include "radiosity.h"
  
  
  
  /***************************************
   *
   *    Global shared variables
   *
   ****************************************/
  
  Global *global;
  Timing **timing;
  
  
  /***************************************
   *
   *    Global variables (not shared)
   *
   ****************************************/
  
  int   n_processors              = DEFAULT_N_PROCESSORS ;
  int   n_taskqueues              = DEFAULT_N_TASKQUEUES ;
  int   n_tasks_per_queue         = DEFAULT_N_TASKS_PER_QUEUE ;
  int   N_inter_parallel_bf_refine= DEFAULT_N_INTER_PARALLEL_BF_REFINEMENT ;
  int   N_visibility_per_task     = DEFAULT_N_VISIBILITY_PER_TASK ;
  float Area_epsilon              = DEFAULT_AREA_EPSILON ;
  float Energy_epsilon            = DEFAULT_ENERGY_EPSILON ;
  float BFepsilon                 = DEFAULT_BFEPSILON ;
  
  int batch_mode = 0 ;
  int verbose_mode = 0 ;
  
  /*
    in converting from a fork process model to an sproc (threads) model,
    taskqueue_id and time_process_start are converted to individual arrays
    without worrying about false sharing.  This is because taskqueue_id is 
    read-only  once written by the parent process, and time_process_start 
    is written only once by each process.
    */
  
  int taskqueue_id[MAX_PROCESSORS] ; 		/* Task queue ID */
  int time_rad_start, time_rad_end, time_process_start[MAX_PROCESSORS] ;
  
  
  /*********************************************************
   *
   *    Global variables (used only by the master process)
   *
   **********************************************************/
  
#define N_SLIDERS (5)
  
  void change_view_x(), change_view_y(), change_view_zoom() ;
  void change_BFepsilon(), change_EXepsilon(), change_area_epsilon() ;
  
  slider sliders[N_SLIDERS] = {
      { "View(X)  deg ", -100,  100, (int)DFLT_VIEW_ROT_X,  5,  change_view_x },
      { "View(Y)  deg ", -100,  100, (int)DFLT_VIEW_ROT_Y,  5,  change_view_y },
      { "View(Zoom)   ",   0,  50, (int)DFLT_VIEW_ZOOM*10,6,  change_view_zoom },
      { "BF-e      0.1%",  0,  50,  (int)(DEFAULT_BFEPSILON *1000.0),
            11, change_BFepsilon },
      { "Area-e       ",   0, 5000, (int)DEFAULT_AREA_EPSILON,
            11, change_area_epsilon },
  } ;
  
#define N_CHOICES (4)
  
#define CHOICE_RAD_RUN    (0)
#define CHOICE_RAD_STEP   (1)
#define CHOICE_RAD_RESET  (2)
  
#define CHOICE_DISP_RADIOSITY   (0)
#define CHOICE_DISP_SHADED      (1)
#define CHOICE_DISP_PATCH       (2)
#define CHOICE_DISP_MESH        (3)
#define CHOICE_DISP_INTERACTION (4)
  
#define CHOICE_UTIL_PS        (0)
#define CHOICE_UTIL_STAT_CRT  (1)
#define CHOICE_UTIL_STAT_FILE (2)
#define CHOICE_UTIL_CLEAR_RAD (3)
  
  void change_display() ;
  void start_radiosity() ;
  void select_model() ;
  void utility_tools() ;
  
  choice choices[N_CHOICES] = {
      { "Run",
            { "Run", "Step", "Reset", 0 },
            0, start_radiosity }, 
      { "Display",
            { "Filled",   "Smooth shading", "Show polygon edges",
                  "Show element edges",  "Show interactions", 0 },
            0, change_display }, 
      { "Models",
            { "Test", "Room", "LargeRoom", 0 },
            0, select_model }, 
      { "Tools",
            { "HardCopy(PS)", "Statistics", "Statistics(file)",
                  "Clear Radiosity Value", 0 },
            0, utility_tools }, 
  } ;
  
  static void expose_callback() ;
  
  /***************************************
   *
   *    Main function.
   *
   ****************************************/
  
  static void parse_args() ;
  
  static int dostats = 0;
  
  main( ac, av )
  
  int ac ;
  char *av[] ;
{
    void radiosity(), parse_args() ;
    int i;
    unsigned long int total_rad_time, max_rad_time, min_rad_time;
    unsigned long int total_refine_time, max_refine_time, min_refine_time;
    unsigned long int total_wait_time, max_wait_time, min_wait_time;
    unsigned long int total_vertex_time, max_vertex_time, min_vertex_time;
    
    /* Parse arguments */
    parse_args( ac, av ) ;
    choices[2].init_value = model_selector ;
    
    /* Initialize graphic device */
    if( batch_mode == 0 )
        {
            g_init( ac, av ) ;
            setup_view( DFLT_VIEW_ROT_X, DFLT_VIEW_ROT_Y,
                       DFLT_VIEW_DIST, DFLT_VIEW_ZOOM,0 ) ;
        }
    
    /* Initialize ANL macro */
    {__tid__[__threads__++]=pthread_self();} ;
    
    /* Allocate global shared memory and initialize */
    global = (Global *) our_malloc(sizeof(Global),__FILE__,__LINE__); ;
    if( global == 0 )
        {
            printf( "Can't allocate memory\n" ) ;
            exit(1) ;
        }
    init_global(0) ;
    
    timing = (Timing **) our_malloc(n_processors * sizeof(Timing *),__FILE__,__LINE__);;
    for (i = 0; i < n_processors; i++)
        timing[i] = (Timing *) our_malloc(sizeof(Timing),__FILE__,__LINE__);;
    
    /* Initialize shared lock */
    init_sharedlock(0) ;
    
    /* Initial random testing rays array for visibility test. */
    init_visibility_module(0) ;
    
/* POSSIBLE ENHANCEMENT:  Here is where one might distribute the 
   sobj_struct, task_struct, and vis_struct data structures across 
   physically distributed memories as desired.

   One way to place data is as follows:

   int i;

   for (i=0;i<n_processors;i++) {
     Place all addresses x such that
       &(sobj_struct[i]) <= x < &(sobj_struct[i+1]) on node i
     Place all addresses x such that
       &(task_struct[i]) <= x < &(task_struct[i+1]) on node i
     Place all addresses x such that
       &(vis_struct[i]) <= x < &(vis_struct[i+1]) on node i
   }

*/

    if( batch_mode )
        {
            /* In batch mode, create child processes and start immediately */
            
            /* Time stamp */
            {long time(); (time_rad_start ) = time(0);};
            
            global->index = 0;
            for( i = 0 ; i < n_processors-1 ; i++ )
                {
                    taskqueue_id[i] = assign_taskq(0) ;
                    {
pthread_mutex_lock(&__intern__);
assert(__threads__<__MAX_THREADS__);
pthread_create(&(__tid__[__threads__++]), NULL, (void*(*)(void *))(radiosity), NULL);
pthread_mutex_unlock(&__intern__);
};
                }
            
            /* And start processing */
            taskqueue_id[n_processors] = assign_taskq(0) ;
            radiosity() ;
            
            {int aantal; for(aantal=n_processors-1;aantal>0;aantal--) pthread_join(__tid__[--__threads__], NULL);};
            
            /* Time stamp */
            {long time(); (time_rad_end ) = time(0);};
            
            /* Print out running time */
            printf("TIMING STATISTICS MEASURED BY MAIN PROCESS:\n");
            
            print_running_time(0);
            
            if (dostats) {
                printf("\n\n\nPER-PROCESS STATISTICS:\n");
                
                printf("%8s%20s%20s%12s%12s\n","Proc","Total","Refine","Wait","Smooth");
                printf("%8s%20s%20s%12s%12s\n\n","","Time","Time","Time","Time");
                for (i = 0; i < n_processors; i++)
                    printf("%8d%20lu%20lu%12lu%12lu\n",i,timing[i]->rad_time, timing[i]->refine_time, timing[i]->wait_time, timing[i]->vertex_time);
                
                total_rad_time = timing[0]->rad_time;
                max_rad_time = timing[0]->rad_time;
                min_rad_time = timing[0]->rad_time;
                
                total_refine_time = timing[0]->refine_time;
                max_refine_time = timing[0]->refine_time;
                min_refine_time = timing[0]->refine_time;
                
                total_wait_time = timing[0]->wait_time;
                max_wait_time = timing[0]->wait_time;
                min_wait_time = timing[0]->wait_time;
                
                total_vertex_time = timing[0]->vertex_time;
                max_vertex_time = timing[0]->vertex_time;
                min_vertex_time = timing[0]->vertex_time;
                
                for (i = 1; i < n_processors; i++) {
                    total_rad_time += timing[i]->rad_time;
                    if (timing[i]->rad_time > max_rad_time)
                        max_rad_time = timing[i]->rad_time;
                    if (timing[i]->rad_time < min_rad_time)
                        min_rad_time = timing[i]->rad_time;
                    
                    total_refine_time += timing[i]->refine_time;
                    if (timing[i]->refine_time > max_refine_time)
                        max_refine_time = timing[i]->refine_time;
                    if (timing[i]->refine_time < min_refine_time)
                        min_refine_time = timing[i]->refine_time;
                    
                    total_wait_time += timing[i]->wait_time;
                    if (timing[i]->wait_time > max_wait_time)
                        max_wait_time = timing[i]->wait_time;
                    if (timing[i]->wait_time < min_wait_time)
                        min_wait_time = timing[i]->wait_time;
                    
                    total_vertex_time += timing[i]->vertex_time;
                    if (timing[i]->vertex_time > max_vertex_time)
                        max_vertex_time = timing[i]->vertex_time;
                    if (timing[i]->vertex_time < min_vertex_time)
                        min_vertex_time = timing[i]->vertex_time;
                }
                
                printf("\n\n%8s%20lu%20lu%12lu%12lu\n","Max", max_rad_time, max_refine_time, max_wait_time, max_vertex_time);
                printf("\n%8s%20lu%20lu%12lu%12lu\n","Min", min_rad_time, min_refine_time, min_wait_time, min_vertex_time);
                printf("\n%8s%20lu%20lu%12lu%12lu\n","Avg", (int) (((double) total_rad_time) / ((double) (1.0 * n_processors))), (int) (((double) total_refine_time) / ((double) (1.0 * n_processors))), (int) (((double) total_wait_time) / ((double) (1.0 * n_processors))), (int) (((double) total_vertex_time) / ((double) (1.0 * n_processors))));
                printf("\n\n");
                
            }
            
            /*	print_fork_time(0) ; */
            
            print_statistics( stdout, 0 ) ;
        }
    else
        {
            /* In interactive mode, start workers, and the master starts
               notification loop */
            
            /* Start notification loop */
            g_start( expose_callback,
                    N_SLIDERS, sliders, N_CHOICES, choices ) ;
        }
    {exit(0);};
    exit(0) ;
}



/***************************************
 *
 *    PANEL call back routine
 *
 *    start_radiosity()   (MASTER only)
 *
 ****************************************/

static int disp_fill_switch = 1 ;
static int disp_shade_switch = 0 ;
static int disp_fill_mode = 1 ;
static int disp_patch_switch = 0 ;
static int disp_mesh_switch  = 0 ;
static int disp_interaction_switch = 0 ;
static int disp_crnt_view_x = (int)DFLT_VIEW_ROT_X ;
  static int disp_crnt_view_y = (int)DFLT_VIEW_ROT_Y ;
  static float disp_crnt_zoom = DFLT_VIEW_ZOOM ;
  
  
#if defined(SGI_GL) && defined(GL_NASA)
  void start_radiosity( ap )
  Actuator *ap ;
#else
  void start_radiosity( val )
  int val ;
#endif
{
    static int state = 0 ;
    void radiosity() ;
    int init_ray_tasks() ;
    void init_radavg_tasks() ;
    int i;
    unsigned total_rad_time, max_rad_time, min_rad_time;
    unsigned total_refine_time, max_refine_time, min_refine_time;
    unsigned total_wait_time, max_wait_time, min_wait_time;
    unsigned total_vertex_time, max_vertex_time, min_vertex_time;
    
#if defined(SGI_GL) && defined(GL_NASA)
    int val ;
    
    val = g_get_choice_val( ap, &choices[0] ) ;
#endif
    
    if( val == CHOICE_RAD_RUN )
        {
            if( state == -1 )
                {
                    printf( "Please reset first\007\n" ) ;
                    return ;
                }
            
            /* Time stamp */
            {long time(); (time_rad_start ) = time(0);} ;
            
            
            global->index = 0;
            
            /* Create slave processes */
            for (i = 0 ; i < n_processors-1 ; i++ )
                {
                    taskqueue_id[i] = assign_taskq(0) ;
                    {
pthread_mutex_lock(&__intern__);
assert(__threads__<__MAX_THREADS__);
pthread_create(&(__tid__[__threads__++]), NULL, (void*(*)(void *))(radiosity), NULL);
pthread_mutex_unlock(&__intern__);
};
                }
            
            /* And start processing */
            taskqueue_id[n_processors] = assign_taskq(0) ;
            radiosity() ;
            
            {int aantal; for(aantal=n_processors-1;aantal>0;aantal--) pthread_join(__tid__[--__threads__], NULL);};
            /* Time stamp */
            {long time(); (time_rad_end ) = time(0);};
            
            /* Print out running time */
            /* Print out running time */
            printf("TIMING STATISTICS MEASURED BY MAIN PROCESS:\n");
            
            print_running_time(0);
            
            if (dostats) {
                printf("\n\n\nPER-PROCESS STATISTICS:\n");
                
                printf("%8s%20s%20s%12s%12s\n","Proc","Total","Refine","Wait","Smooth");
                printf("%8s%20s%20s%12s%12s\n\n","","Time","Time","Time","Time")
                    ;
                for (i = 0; i < n_processors; i++)
                    printf("%8d%20lu%20lu%12lu%12lu\n",i,timing[i]->rad_time, timing[i]->refine_time, timing[i]->wait_time, timing[i]->vertex_time);
                
                total_rad_time = timing[0]->rad_time;
                max_rad_time = timing[0]->rad_time;
                min_rad_time = timing[0]->rad_time;
                
                total_refine_time = timing[0]->refine_time;
                max_refine_time = timing[0]->refine_time;
                min_refine_time = timing[0]->refine_time;
                
                total_wait_time = timing[0]->wait_time;
                max_wait_time = timing[0]->wait_time;
                min_wait_time = timing[0]->wait_time;
                
                total_vertex_time = timing[0]->vertex_time;
                max_vertex_time = timing[0]->vertex_time;
                min_vertex_time = timing[0]->vertex_time;
                
                for (i = 1; i < n_processors; i++) {
                    total_rad_time += timing[i]->rad_time;
                    if (timing[i]->rad_time > max_rad_time)
                        max_rad_time = timing[i]->rad_time;
                    if (timing[i]->rad_time < min_rad_time)
                        min_rad_time = timing[i]->rad_time;
                    
                    total_refine_time += timing[i]->refine_time;
                    if (timing[i]->refine_time > max_refine_time)
                        max_refine_time = timing[i]->refine_time;
                    if (timing[i]->refine_time < min_refine_time)
                        min_refine_time = timing[i]->refine_time;
                    
                    total_wait_time += timing[i]->wait_time;
                    if (timing[i]->wait_time > max_wait_time)
                        max_wait_time = timing[i]->wait_time;
                    if (timing[i]->wait_time < min_wait_time)
                        min_wait_time = timing[i]->wait_time;
                    
                    total_vertex_time += timing[i]->vertex_time;
                    if (timing[i]->vertex_time > max_vertex_time)
                        max_vertex_time = timing[i]->vertex_time;
                    if (timing[i]->vertex_time < min_vertex_time)
                        min_vertex_time = timing[i]->vertex_time;
                }
                
                
                printf("\n\n%8s%20lu%20lu%12lu%12lu\n","Max", max_rad_time, max_refine_time, max_wait_time, max_vertex_time);
                printf("\n%8s%20lu%20lu%12lu%12lu\n","Min", min_rad_time, min_refine_time, min_wait_time, min_vertex_time);
                printf("\n%8s%20lu%20lu%12lu%12lu\n","Avg", (int) (((double) total_rad_time) / ((double) (1.0 * n_processors))), (int) (((double) total_refine_time) / ((double) (1.0 * n_processors))), (int) (((double) total_wait_time) / ((double) (1.0 * n_processors))), (int) (((double) total_vertex_time) / ((double) (1.0 * n_processors))));
                printf("\n\n");
                
            }
            
            /*      print_fork_time(0) ; */
            
            print_statistics( stdout, 0 ) ;
            
            /* Display image */
            display_scene( disp_fill_mode, disp_patch_switch,
                          disp_mesh_switch, disp_interaction_switch, 0) ;
            
            state = -1 ;
        }
    
    else if( val == CHOICE_RAD_STEP )
        {
            if( state == -1 )
                {
                    printf( "Please reset first\007\n" ) ;
                    return ;
                }
            
            /* Step execution */
            switch( state )
                {
                case 0:
                    /* Step execute as a single process */
                    
                    global->index = 1;
                    /* Create slave processes */
                    for ( i = 0 ; i < n_processors-1 ; i++ )
                        {
                            taskqueue_id[i] = assign_taskq(0) ;
                            {
pthread_mutex_lock(&__intern__);
assert(__threads__<__MAX_THREADS__);
pthread_create(&(__tid__[__threads__++]), NULL, (void*(*)(void *))(radiosity), NULL);
pthread_mutex_unlock(&__intern__);
};
                        }
                    
                    taskqueue_id[n_processors] = assign_taskq(0) ;
                    
                    /* Decompose model objects into patches and build
                       the BSP tree */
                    /* Create the first tasks (MASTER only) */
                    init_modeling_tasks(0) ;
                    process_tasks(0) ;
                    state ++ ;
                    break ;
                    
                case 1:
                    if( init_ray_tasks(0) )
                        {
                            {
 pthread_mutex_lock(&((global->barrier).bar_mutex));
 (global->barrier).bar_teller++;
 if ((global->barrier).bar_teller == (n_processors)) {
     (global->barrier).bar_teller = 0;
     pthread_cond_broadcast(&((global->barrier).bar_cond));
 } else
     pthread_cond_wait(&((global->barrier).bar_cond), &((global->barrier).bar_mutex));
 pthread_mutex_unlock(&((global->barrier).bar_mutex));
};
                            process_tasks(0) ;
                        }
                    else
                        state++ ;
                    break ;
                default:
                    {
 pthread_mutex_lock(&((global->barrier).bar_mutex));
 (global->barrier).bar_teller++;
 if ((global->barrier).bar_teller == (n_processors)) {
     (global->barrier).bar_teller = 0;
     pthread_cond_broadcast(&((global->barrier).bar_cond));
 } else
     pthread_cond_wait(&((global->barrier).bar_cond), &((global->barrier).bar_mutex));
 pthread_mutex_unlock(&((global->barrier).bar_mutex));
};
                    init_radavg_tasks( RAD_AVERAGING_MODE, 0 ) ;
                    process_tasks(0) ;
                    init_radavg_tasks( RAD_NORMALIZING_MODE, 0 ) ;
                    process_tasks(0) ;
                    
                    {int aantal; for(aantal=n_processors-1;aantal>0;aantal--) pthread_join(__tid__[--__threads__], NULL);}
                        state = -1 ;
                }
            
            /* Display image */
            display_scene( disp_fill_mode, disp_patch_switch,
                          disp_mesh_switch, disp_interaction_switch, 0) ;
        }
    
    else if( val == CHOICE_RAD_RESET )
        {
            /* Initialize global variables again */
            init_global(0) ;
            init_visibility_module(0) ;
            g_clear() ;
            state = 0 ;
        }
}


/***************************************
 *
 *    PANEL call back routine
 *
 *    change_display()   (MASTER only)
 *
 ****************************************/


#if defined(SGI_GL) && defined(GL_NASA)
void change_display( ap )
  
  Actuator *ap ;
#else
  void change_display( val )
  
  int val ;
#endif
{
#if defined(SGI_GL) && defined(GL_NASA)
    int val ;
    
    val = g_get_choice_val( ap, &choices[1] ) ;
#endif
    
    /* Display image */
    switch( val )
        {
        case CHOICE_DISP_RADIOSITY:
            disp_fill_switch = (! disp_fill_switch) ;
            break ;
        case CHOICE_DISP_SHADED:
            disp_shade_switch = (! disp_shade_switch) ;
            break ;
        case CHOICE_DISP_PATCH:
            disp_patch_switch = (! disp_patch_switch) ;
            break ;
        case CHOICE_DISP_MESH:
            disp_mesh_switch = (! disp_mesh_switch) ;
            break ;
        case CHOICE_DISP_INTERACTION:
            disp_interaction_switch = (! disp_interaction_switch) ;
            break ;
        default:
            return ;
        }
    
    if( disp_fill_switch == 0 )
        disp_fill_mode = 0 ;
    else
        {
            if( disp_shade_switch == 0 )
                disp_fill_mode = 1 ;
            else
                disp_fill_mode = 2 ;
        }
    
    /* Display image */
    display_scene( disp_fill_mode, disp_patch_switch,
                  disp_mesh_switch, disp_interaction_switch, 0 ) ;
}


/*****************************************************
 *
 *    PANEL call back routine
 *
 *    change_view_y()            (MASTER only)
 *    change_BFepsilon()
 *    change_area_epsilon()
 *
 ******************************************************/

static void change_view()
{
    /* Change the view */
    setup_view( (float)disp_crnt_view_x, (float)disp_crnt_view_y,
               DFLT_VIEW_DIST, disp_crnt_zoom, 0 ) ;
    
    /* And redraw */
    display_scene( disp_fill_mode, disp_patch_switch,
                  disp_mesh_switch, disp_interaction_switch, 0 ) ;
}


#if defined(SGI_GL) && defined(GL_NASA)
void change_view_x( ap )
  Actuator *ap ;
#else
  void change_view_x( val )
  int val ;
#endif
{
#if defined(SGI_GL) && defined(GL_NASA)
    int val = g_get_slide_val( ap ) ;
#endif
    
    /* Save current rot-X value */
    disp_crnt_view_x = val ;
    change_view() ;
}


#if defined(SGI_GL) && defined(GL_NASA)
void change_view_y( ap )
  Actuator *ap ;
#else
  void change_view_y( val )
  int val ;
#endif
{
#if defined(SGI_GL) && defined(GL_NASA)
    int val = g_get_slide_val( ap ) ;
#endif
    
    /* Save current rot-Y value */
    disp_crnt_view_y = val ;
    change_view() ;
}


#if defined(SGI_GL) && defined(GL_NASA)
void change_view_zoom( ap )
  Actuator *ap ;
#else
  void change_view_zoom( val )
  int val ;
#endif
{
#if defined(SGI_GL) && defined(GL_NASA)
    int val = g_get_slide_val( ap ) ;
#endif
    
    /* Save current zoom value */
    disp_crnt_zoom = (float)val / 10.0 ;
    change_view() ;
}


#if defined(SGI_GL) && defined(GL_NASA)
void change_BFepsilon( ap )
  Actuator *ap ;
#else
  void change_BFepsilon( val )
  int val ;
#endif
{
#if defined(SGI_GL) && defined(GL_NASA)
    int val = g_get_slide_val( ap ) ;
#endif
    BFepsilon = (float)val / 1000.0 ;
}



#if defined(SGI_GL) && defined(GL_NASA)
void change_area_epsilon( ap )
  Actuator *ap ;
#else
  void change_area_epsilon( val )
  int val ;
#endif
{
#if defined(SGI_GL) && defined(GL_NASA)
    int val = g_get_slide_val( ap ) ;
#endif
    Area_epsilon = (float)val ;
}


/***************************************
 *
 *    PANEL call back routine
 *
 *    select_model()   (MASTER only)
 *
 ****************************************/

#if defined(SGI_GL) && defined(GL_NASA)
void select_model( ap )
  Actuator *ap ;
#else
  void select_model( val )
  int val ;
#endif
{
#if defined(SGI_GL) && defined(GL_NASA)
    int val = g_get_choice_val( ap, &choices[2] ) ;
#endif
    switch( val )
        {
        case MODEL_TEST_DATA:
            model_selector = MODEL_TEST_DATA ;
            break ;
        case MODEL_ROOM_DATA:
            model_selector = MODEL_ROOM_DATA ;
            break ;
        case MODEL_LARGEROOM_DATA:
            model_selector = MODEL_LARGEROOM_DATA ;
            break ;
        }
}



/***************************************
 *
 *    PANEL call back routine
 *
 *    utility_tools()   (MASTER only)
 *
 ****************************************/


#if defined(SGI_GL) && defined(GL_NASA)
void utility_tools( ap )
  Actuator *ap ;
#else
  void utility_tools( val )
  int val ;
#endif
{
    FILE *fd ;
#if defined(SGI_GL) && defined(GL_NASA)
    int val = g_get_choice_val( ap, &choices[3] ) ;
#endif
    
    switch( val )
        {
        case CHOICE_UTIL_PS:
            /* Open PS file */
            ps_open( "radiosity.ps" ) ;
            
            /* Change the view */
            ps_setup_view( DFLT_VIEW_ROT_X, (float)disp_crnt_view_y,
                          DFLT_VIEW_DIST, DFLT_VIEW_ZOOM, 0 ) ;
            
            /* And redraw */
            ps_display_scene( disp_fill_mode, disp_patch_switch,
                             disp_mesh_switch, disp_interaction_switch, 0 ) ;
            
            /* Close */
            ps_close() ;
            break ;
        case CHOICE_UTIL_STAT_CRT:
            print_statistics( stdout, 0 ) ;
            break ;
        case CHOICE_UTIL_STAT_FILE:
            if( (fd = fopen( "radiosity_stat", "w" )) == 0 )
                {
                    perror( "radiosity_stat" ) ;
                    break ;
                }
            print_statistics( fd, 0 ) ;
            fclose( fd ) ;
            break ;
        case CHOICE_UTIL_CLEAR_RAD:
            clear_radiosity(0) ;
        }
}


/***************************************
 *
 *    Exposure call back
 *
 ****************************************/

static void expose_callback()
{
    /* Display image */
    display_scene( disp_fill_mode, disp_patch_switch,
                  disp_mesh_switch, disp_interaction_switch, 0 ) ;
}


/***************************************
 *
 *    radiosity()  Radiosity task main
 *
 ****************************************/


void radiosity()
{
    int init_ray_tasks() ;
    unsigned process_id;
    void init_radavg_tasks() ;
    unsigned long int rad_start, refine_done, vertex_start, vertex_done;
    
    {pthread_mutex_lock(&(global->index_lock));};
    process_id = global->index++;
    {pthread_mutex_unlock(&(global->index_lock));};
    process_id = process_id % n_processors;
    
    if ((process_id == 0) || (dostats))
        {long time(); (rad_start) = time(0);};
    
    /* POSSIBLE ENHANCEMENT:  Here is where one might pin processes to
       processors to avoid migration */

    /* POSSIBLE ENHANCEMENT:  Here is where one might reset the
       statistics that one is measuring about the parallel execution */

    /* Decompose model objects into patches and build the BSP tree */
    /* Create the initial tasks */
    init_modeling_tasks(process_id) ;
    process_tasks(process_id) ;
    
    /* Gather rays & do BF refinement */
    while( init_ray_tasks(process_id) )
        {
            /* Wait till tasks are put in the queue */
            {
 pthread_mutex_lock(&((global->barrier).bar_mutex));
 (global->barrier).bar_teller++;
 if ((global->barrier).bar_teller == (n_processors)) {
     (global->barrier).bar_teller = 0;
     pthread_cond_broadcast(&((global->barrier).bar_cond));
 } else
     pthread_cond_wait(&((global->barrier).bar_cond), &((global->barrier).bar_mutex));
 pthread_mutex_unlock(&((global->barrier).bar_mutex));
};
            /* Then perform ray-gathering and BF-refinement till the
               solution converges */
            process_tasks(process_id) ;
        }
    
    if ((process_id == 0) || (dostats))
        {long time(); (refine_done) = time(0);};
    
    {
 pthread_mutex_lock(&((global->barrier).bar_mutex));
 (global->barrier).bar_teller++;
 if ((global->barrier).bar_teller == (n_processors)) {
     (global->barrier).bar_teller = 0;
     pthread_cond_broadcast(&((global->barrier).bar_cond));
 } else
     pthread_cond_wait(&((global->barrier).bar_cond), &((global->barrier).bar_mutex));
 pthread_mutex_unlock(&((global->barrier).bar_mutex));
};
    
    if ((process_id == 0) || (dostats))
        {long time(); (vertex_start) = time(0);};
    
    /* Compute area-weighted radiosity value at each vertex */
    init_radavg_tasks( RAD_AVERAGING_MODE, process_id ) ;
    process_tasks(process_id) ;
    
    /* Then normalize the radiosity at vertices */
    init_radavg_tasks( RAD_NORMALIZING_MODE, process_id ) ;
    process_tasks(process_id) ;
    
    if ((process_id == 0) || (dostats))
        {long time(); (vertex_done) = time(0);};
    
    if ((process_id == 0) || (dostats)) {
        timing[process_id]->rad_start = rad_start;
        timing[process_id]->rad_time = vertex_done - rad_start;
        timing[process_id]->refine_time = refine_done - rad_start;
        timing[process_id]->vertex_time = vertex_done - vertex_start;
        timing[process_id]->wait_time = vertex_start - refine_done;
    }
    
}



/***************************************************************************
 *
 *    init_ray_tasks()
 *
 *    Create initial tasks to perform ray gathering.
 *
 ****************************************************************************/


#if PATCH_ASSIGNMENT == PATCH_ASSIGNMENT_STATIC
static void _init_ray_tasks_static() ;
#define _INIT_RAY_TASK  _init_ray_tasks_static
#endif
  
  
  static void _init_ray_tasks_cost2() ;
  static int avg_cost_of_q ;
  static int avg_cost_of_patch ;
  static int cost_of_this_q ;
  static int crnt_qid ;
  static int queue_cost[MAX_TASKQUEUES] ;
  
#if PATCH_ASSIGNMENT == PATCH_ASSIGNMENT_COSTBASED
#define _INIT_RAY_TASK  _init_ray_tasks_cost2
#endif
  
  
  
  
  int init_ray_tasks(process_id)
  unsigned process_id;
{
    int conv ;
    
    /* If this is not the first process to initialize, then return */
    {pthread_mutex_lock(&(global->avg_radiosity_lock));};
    if( ! check_task_counter(process_id) )
        {
            conv = global->converged ;
            {pthread_mutex_unlock(&(global->avg_radiosity_lock));};
            return( conv == 0 ) ;
        }
    
    /* Check radiosity convergence */
    conv = radiosity_converged(process_id) ;
    global->converged = conv ;
    
    /* Reset total energy variable */
    global->prev_total_energy = global->total_energy ;
    global->total_energy.r = 0.0 ;
    global->total_energy.g = 0.0 ;
    global->total_energy.b = 0.0 ;
    global->total_patch_area = 0.0 ;
    
    /* Increment iteration counter */
    global->iteration_count++ ;
    {pthread_mutex_unlock(&(global->avg_radiosity_lock));};
    
    /* If radiosity converged, then return 0 */
    if( conv )
        return( 0 ) ;
    
    
#if PATCH_ASSIGNMENT == PATCH_ASSIGNMENT_COSTBASED
    /* Compute average cost per queue. Also reset the cost variable.
       The 'cost_sum' is not locked since no one is processing rays
       at this moment */
    
    for( crnt_qid = 0 ; crnt_qid < n_taskqueues ; crnt_qid++ )
        queue_cost[ crnt_qid ] = 0 ;
    
    avg_cost_of_q = global->cost_estimate_sum / n_taskqueues ;
    avg_cost_of_patch = global->cost_estimate_sum / global->n_total_patches ;
    cost_of_this_q = 0 ;
    crnt_qid = 0 ;
    
    global->cost_sum = 0 ;
    global->cost_estimate_sum = 0 ;
    
    /* layered selection of tasks */
    foreach_patch_in_bsp( _INIT_RAY_TASK, 2, process_id ) ;
    foreach_patch_in_bsp( _INIT_RAY_TASK, 1, process_id ) ;
#endif
    
    /* Create BF refinement tasks */
    foreach_patch_in_bsp( _INIT_RAY_TASK, 0, process_id ) ;
    
    return( 1 ) ;
}


static void _init_ray_tasks_static( p, dummy, process_id )
  Patch *p ;
  int dummy ;
  unsigned process_id;
{
    /* Clear incoming energy variable */
    p->el_root->rad_in.r = 0.0 ;
    p->el_root->rad_in.g = 0.0 ;
    p->el_root->rad_in.b = 0.0 ;
    
    enqueue_ray_task( (p->seq_no >> 2) % n_taskqueues, p->el_root,
                     TASK_APPEND, process_id ) ;
}


static void _init_ray_tasks_cost2( p, layer, process_id )
  Patch *p ;
  int layer ;
  unsigned process_id;
{
    Patch_Cost *pc ;
    int c_est ;
    int qid ;
    int min_cost_q, min_cost ;
    
    
    pc = &global->patch_cost[ p->seq_no ] ;
    c_est = pc->cost_estimate ;
    
    if( c_est < 0 )
        /* Already processed */
        return ;
    
    if( c_est < avg_cost_of_patch * layer )
        return ;
    
    /* Find the first available queue */
    min_cost_q = 0 ;
    min_cost = queue_cost[ 0 ] ;
    for( qid = 0 ; qid < n_taskqueues ; qid++ )
        {
            if( (c_est + queue_cost[ qid ]) <= avg_cost_of_q )
                break ;
            
            if( min_cost > queue_cost[ qid ] )
                {
                    min_cost_q = qid ;
                    min_cost = queue_cost[ qid ] ;
                }
        }
    
    if( qid >= n_taskqueues )
        {
            /* All queues are nearly full. Put to min-cost queue */
            qid = min_cost_q ;
        }
    
    /* Update queue cost */
    queue_cost[ qid ] += c_est ;
    
    
    /* Clear incoming energy variable */
    p->el_root->rad_in.r = 0.0 ;
    p->el_root->rad_in.g = 0.0 ;
    p->el_root->rad_in.b = 0.0 ;
    
    /* Clear cost value */
    pc->cost_estimate = -1 ;
    pc->n_bsp_node    = 0 ;
    
    /* Enqueue */
    enqueue_ray_task( qid, p->el_root, TASK_APPEND, process_id ) ;
    
}




/***************************************************************************
 *
 *    init_radavg_tasks()
 *
 *    Create initial tasks to perform radiosity averaging.
 *
 ****************************************************************************/

static void _init_radavg_tasks() ;
  
  void init_radavg_tasks( mode, process_id )
  
  int mode ;
  unsigned process_id;
{
    
    /* If this is not the first process to initialize, then return */
    if( ! check_task_counter(process_id ) )
        return ;
    
    /* Create RadAvg tasks */
    foreach_patch_in_bsp( _init_radavg_tasks, mode, process_id ) ;
}


static void _init_radavg_tasks( p, mode, process_id )
  Patch *p ;
  int mode ;
  unsigned process_id;
{
    enqueue_radavg_task( p->seq_no % n_taskqueues, p->el_root, mode, process_id  ) ;
}



/***************************************
 *
 *    init_global()
 *
 ****************************************/


void init_global(process_id)
  unsigned process_id;
{
    /* Clear BSP root pointer */
    global->index = 1;  /* ****** */
    global->bsp_root = 0 ;
    {pthread_mutex_init(&(global->index_lock),NULL);};
    {pthread_mutex_init(&(global->bsp_tree_lock),NULL);};
    
    /* Initialize radiosity statistics variables */
    {pthread_mutex_init(&(global->avg_radiosity_lock),NULL);};
    global->converged = 0 ;
    global->prev_total_energy.r = 0.0 ;   
    global->prev_total_energy.g = 0.0 ;   
    global->prev_total_energy.b = 0.0 ;   
    global->total_energy.r = 1.0 ;   
    global->total_energy.g = 1.0 ;   
    global->total_energy.b = 1.0 ;   
    global->total_patch_area = 1.0 ;
    global->iteration_count = -1 ;     /* init_ray_task() increments to 0 */
    
    /* Initialize the cost sum */
    {pthread_mutex_init(&(global->cost_sum_lock),NULL);};
    global->cost_sum = 0 ;
    global->cost_estimate_sum = 0 ;
    
    /* Initialize the barrier */
    {
pthread_mutex_init(&((global->barrier).bar_mutex), NULL);
pthread_cond_init(&((global->barrier).bar_cond), NULL);
(global->barrier).bar_teller=0;
};
    {pthread_mutex_init(&(global->pbar_lock),NULL);};
    global->pbar_count = 0 ;
    
    /* Initialize task counter */
    global->task_counter = 0 ;
    {pthread_mutex_init(&(global->task_counter_lock),NULL);};
    
    /* Initialize task queue */
    init_taskq(process_id) ;
    
    /* Initialize Patch, Element, Interaction free lists */
    init_patchlist(process_id) ;
    init_elemlist(process_id) ;
    init_interactionlist(process_id) ;
    init_elemvertex(process_id) ;
    init_edge(process_id) ;
    
    /* Initialize statistical info */
    init_stat_info(process_id) ;
    
}


/*************************************************************
 *
 * parse_args()   Parse arguments
 *
 **************************************************************/

static void parse_args( ac, av )
  
  int ac ;
  char *av[] ;
{
    int cnt ;
    
    /* Parse arguments */
    for( cnt = 1 ; cnt < ac ; cnt++ ) 
        {
            if( strcmp( av[cnt], "-p" ) == 0 ) {
                sscanf( av[++cnt], "%d", &n_processors ) ;
                n_taskqueues = n_processors;
            }	
            else if( strcmp( av[cnt], "-tq" ) == 0 ) 
                sscanf( av[++cnt], "%d", &n_tasks_per_queue ) ;
            else if( strcmp( av[cnt], "-ae" ) == 0 )
                sscanf( av[++cnt], "%f", &Area_epsilon ) ;
            else if( strcmp( av[cnt], "-pr" ) == 0 ) 
                sscanf( av[++cnt], "%d", &N_inter_parallel_bf_refine ) ;
            else if( strcmp( av[cnt], "-pv" ) == 0 ) 
                sscanf( av[++cnt], "%d", &N_visibility_per_task ) ;
            else if( strcmp( av[cnt], "-bf" ) == 0 )
                sscanf( av[++cnt], "%f", &BFepsilon ) ;
            else if( strcmp( av[cnt], "-en" ) == 0 )
                sscanf( av[++cnt], "%f", &Energy_epsilon ) ;
            
            else if( strcmp( av[cnt], "-batch" ) == 0 )
                batch_mode = 1 ;
            else if( strcmp( av[cnt], "-verbose" ) == 0 )
                verbose_mode = 1 ;
            else if( strcmp( av[cnt], "-s" ) == 0 )
                dostats = 1 ;
            else if( strcmp( av[cnt], "-room" ) == 0 )
                model_selector = MODEL_ROOM_DATA ;
            else if( strcmp( av[cnt], "-largeroom" ) == 0 )
                model_selector = MODEL_LARGEROOM_DATA ;
            else if(( strcmp( av[cnt], "-help" ) == 0 ) || ( strcmp( av[cnt], "-h" ) == 0 ) || ( strcmp( av[cnt], "-H" ) == 0 ))	    {
                print_usage() ;
                exit(0) ;
            }
        }
    
    
    /* Then check the arguments */
    if( (n_processors < 1) || (MAX_PROCESSORS < n_processors) )
        {
            fprintf( stderr, "Bad number of processors: %d\n",
                    n_processors ) ;
            exit(1) ;
        }
    if( (n_taskqueues < 1) || (MAX_TASKQUEUES < n_taskqueues) )
        {
            fprintf( stderr, "Bad number of task queues: %d\n",
                    n_taskqueues ) ;
            exit(1) ;
        }
    /* Check epsilon values */
    if( Area_epsilon < 0.0 )
        {
            fprintf( stderr, "Area epsilon must be positive\n" ) ;
            exit(1) ;
        }
    if( BFepsilon < 0.0 )
        {
            fprintf( stderr, "BFepsilon must be within [0,1]\n" ) ;
            exit(1) ;
        }
}



/*************************************************************
 *
 *   print_usage()
 *
 **************************************************************/

void print_usage()
{
    fprintf( stderr, "Usage:  RADIOSITY  [options..]\n\n" ) ;
    fprintf( stderr, "\tNote: Must have a space between option label and numeric value, if any\n\n");
    fprintf( stderr, "   -p    (d)  # of processes\n" ) ;
    fprintf( stderr, "   -tq   (d)  # of tasks per queue: default (200) in code for SPLASH\n" ) ;
    fprintf( stderr, "   -ae   (f)  Area epsilon: default (2000.0) in code for SPLASH\n" ) ;
    fprintf( stderr, "   -pr   (d)  # of inter for parallel refinement: default (5) in code for SPLASH\n") ;
    fprintf( stderr, "   -pv   (d)  # of visibility comp in a task: default (4) in code for SPLASH\n") ;
    fprintf( stderr, "   -bf   (f)  BFepsilon (BF refinement): default (0.015) in code for SPLASH\n" ) ;
    fprintf( stderr, "   -en   (f)  Energy epsilon (convergence): default (0.005) in code for SPLASH\n" ) ;
    fprintf( stderr, "   -room      Use room model (default=test)\n" ) ;
    fprintf( stderr, "   -largeroom Use large room model\n" ) ;
    fprintf( stderr, "   -batch     Batch mode (use for SPLASH)\n" ) ;
    fprintf( stderr, "   -verbose   Verbose mode (don't use for SPLASH)\n" ) ;
    fprintf( stderr, "   -s   Measure per-process timing (don't use for SPLASH)\n" ) ;
}


/* Generated from ../Changed/rad_main.C */
