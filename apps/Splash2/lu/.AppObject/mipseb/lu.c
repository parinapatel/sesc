

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

/*************************************************************************/
/*                                                                       */
/*  Parallel dense blocked LU factorization (no pivoting)                */
/*                                                                       */
/*  This version contains two dimensional arrays in which the first      */
/*  dimension is the block to be operated on, and the second contains    */
/*  all data points in that block.  In this manner, all data points in   */
/*  a block (which are operated on by the same processor) are allocated  */
/*  contiguously and locally, and false sharing is eliminated.           */
/*                                                                       */
/*  Command line options:                                                */
/*                                                                       */
/*  -nN : Decompose NxN matrix.                                          */
/*  -pP : P = number of processors.                                      */
/*  -bB : Use a block size of B. BxB elements should fit in cache for    */
/*        good performance. Small block sizes (B=8, B=16) work well.     */
/*  -s  : Print individual processor timing statistics.                  */
/*  -t  : Test output.                                                   */
/*  -o  : Print out matrix values.                                       */
/*  -h  : Print out command line options.                                */
/*                                                                       */
/*  Note: This version works under both the FORK and SPROC models        */
/*                                                                       */
/*************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#if (defined NO_PADDING)
#define PAGE_SIZE                  8
#endif


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


#define MAXRAND                         32767.0
#define DEFAULT_N                         512
#define DEFAULT_P                           1
#if (defined NO_PADDING)
#define DEFAULT_B                           1
#else
#define DEFAULT_B                          16
#endif
#define min(a,b) ((a) < (b) ? (a) : (b))

struct GlobalMemory {
  double *t_in_fac;   
  double *t_in_solve;
  double *t_in_mod; 
  double *t_in_bar;
  double *completion;
  unsigned int starttime; 
  unsigned int rf; 
  unsigned int rs; 
  unsigned int done;
  int id;
  
struct { pthread_mutex_t bar_mutex; pthread_cond_t bar_cond; unsigned bar_teller; } start;

  pthread_mutex_t idlock;
} *Global;

struct LocalCopies {
  double t_in_fac;
  double t_in_solve;
  double t_in_mod;
  double t_in_bar;
};

int n = DEFAULT_N;          /* The size of the matrix */
int P = DEFAULT_P;          /* Number of processors */
int block_size = DEFAULT_B; /* Block dimension */
int nblocks;                /* Number of blocks in each dimension */
int num_rows;               /* Number of processors per row of processor grid */
int num_cols;               /* Number of processors per col of processor grid */
double **a;                 /* a = lu; l and u both placed back in a */
double *rhs;
int *proc_bytes;            /* Bytes to malloc per processor to hold blocks 
			       of A*/
double **last_malloc;       /* Starting point of last block of A */

int test_result = 0;        /* Test result of factorization? */
int doprint = 0;            /* Print out matrix values? */
int dostats = 0;            /* Print out individual processor statistics? */

void SlaveStart();
void OneSolve(int, int, double **, int, int);
void lu0(double *,int, int);
void bdiv(double *, double *, int, int, int, int);
void bmodd(double *, double*, int, int, int, int);
void bmod(double *, double *, double *, int, int, int, int, int, int);
void daxpy(double *, double *, int, double);
int BlockOwner(int, int);
void lu(int, int, int, struct LocalCopies *, int);
void InitA(double *);
double TouchA(int, int);
void PrintA();
void CheckResult(int, double **, double *);
void printerr(char *);


main(argc, argv)

int argc;
char *argv[];

{
  int i, j;
  int ch;
  extern char *optarg;
  int MyNum=0;
  double mint, maxt, avgt;
  double min_fac, min_solve, min_mod, min_bar;
  double max_fac, max_solve, max_mod, max_bar;
  double avg_fac, avg_solve, avg_mod, avg_bar;
  int last_page;
  int proc_num;
  int edge;
  int size;
  unsigned int start;

  {long time(); (start) = time(0);}

  while ((ch = getopt(argc, argv, "n:p:b:cstoh")) != -1) {
    switch(ch) {
    case 'n': n = atoi(optarg); break;
    case 'p': P = atoi(optarg); break;
    case 'b': block_size = atoi(optarg); break;
    case 's': dostats = 1; break;
    case 't': test_result = !test_result; break;
    case 'o': doprint = !doprint; break;
    case 'h': printf("Usage: LU <options>\n\n");
	      printf("options:\n");
              printf("  -nN : Decompose NxN matrix.\n");
              printf("  -pP : P = number of processors.\n");
              printf("  -bB : Use a block size of B. BxB elements should fit in cache for \n");
              printf("        good performance. Small block sizes (B=8, B=16) work well.\n");
              printf("  -c  : Copy non-locally allocated blocks to local memory before use.\n");
              printf("  -s  : Print individual processor timing statistics.\n");
              printf("  -t  : Test output.\n");
              printf("  -o  : Print out matrix values.\n");
              printf("  -h  : Print out command line options.\n\n");
              printf("Default: LU -n%1d -p%1d -b%1d\n",
		     DEFAULT_N,DEFAULT_P,DEFAULT_B);
              exit(0);
              break;
    }
  }

  {__tid__[__threads__++]=pthread_self();}

  printf("\n");
  printf("Blocked Dense LU Factorization\n");
  printf("     %d by %d Matrix\n",n,n);
  printf("     %d Processors\n",P);
  printf("     %d by %d Element Blocks\n",block_size,block_size);
  printf("\n");
  printf("\n");

  num_rows = (int) sqrt((double) P);
  for (;;) {
    num_cols = P/num_rows;
    if (num_rows*num_cols == P)
      break;
    num_rows--;
  }
  nblocks = n/block_size;
  if (block_size * nblocks != n) {
    nblocks++;
  }

  edge = n%block_size;
  if (edge == 0) {
    edge = block_size;
  }
  proc_bytes = (int *) malloc(P*sizeof(int));
  last_malloc = (double **) our_malloc(P*sizeof(double *),__FILE__,__LINE__);;
  for (i=0;i<P;i++) {
    proc_bytes[i] = 0;
    last_malloc[i] = NULL;
  }
  for (i=0;i<nblocks;i++) {
    for (j=0;j<nblocks;j++) {
      proc_num = BlockOwner(i,j);
      if ((i == nblocks-1) && (j == nblocks-1)) {
        size = edge*edge;
      } else if ((i == nblocks-1) || (j == nblocks-1)) {
        size = edge*block_size;
      } else {
        size = block_size*block_size;
      }
      proc_bytes[proc_num] += size*sizeof(double);
    }
  }
  for (i=0;i<P;i++) {
    last_malloc[i] = (double *) our_malloc(proc_bytes[i] + PAGE_SIZE,__FILE__,__LINE__);
    if (last_malloc[i] == NULL) {
      fprintf(stderr,"Could not malloc memory blocks for proc %d\n",i);
      exit(-1);
    } 
/*
    last_malloc[i] = (double *) (((unsigned) last_malloc[i]) + PAGE_SIZE -
                     ((unsigned) last_malloc[i]) % PAGE_SIZE);
*/
    last_malloc[i] = (double *) (((uintptr_t) last_malloc[i]) + PAGE_SIZE -
                     ((uintptr_t) last_malloc[i]) % PAGE_SIZE);

/* Note that this causes all blocks to start out page-aligned, and that
   for block sizes that exceed cache line size, blocks start at cache-line
   aligned addresses as well.  This reduces false sharing */

  }
  a = (double **) our_malloc(nblocks*nblocks*sizeof(double *),__FILE__,__LINE__);;
  if (a == NULL) {
    printerr("Could not malloc memory for a\n");
    exit(-1);
  } 
  for (i=0;i<nblocks;i++) {
    for (j=0;j<nblocks;j++) {
      proc_num = BlockOwner(i,j);
      a[i+j*nblocks] = last_malloc[proc_num];
      if ((i == nblocks-1) && (j == nblocks-1)) {
        size = edge*edge;
      } else if ((i == nblocks-1) || (j == nblocks-1)) {
        size = edge*block_size;
      } else {
        size = block_size*block_size;
      }
      last_malloc[proc_num] += size;
    }
  }

  rhs = (double *) our_malloc(n*sizeof(double),__FILE__,__LINE__);;
  if (rhs == NULL) {
    printerr("Could not malloc memory for rhs\n");
    exit(-1);
  } 

  Global = (struct GlobalMemory *) our_malloc(sizeof(struct GlobalMemory),__FILE__,__LINE__);;
  Global->t_in_fac = (double *) our_malloc(P*sizeof(double),__FILE__,__LINE__);;
  Global->t_in_mod = (double *) our_malloc(P*sizeof(double),__FILE__,__LINE__);;
  Global->t_in_solve = (double *) our_malloc(P*sizeof(double),__FILE__,__LINE__);;
  Global->t_in_bar = (double *) our_malloc(P*sizeof(double),__FILE__,__LINE__);;
  Global->completion = (double *) our_malloc(P*sizeof(double),__FILE__,__LINE__);;

  if (Global == NULL) {
    printerr("Could not malloc memory for Global\n");
    exit(-1);
  } else if (Global->t_in_fac == NULL) {
    printerr("Could not malloc memory for Global->t_in_fac\n");
    exit(-1);
  } else if (Global->t_in_mod == NULL) {
    printerr("Could not malloc memory for Global->t_in_mod\n");
    exit(-1);
  } else if (Global->t_in_solve == NULL) {
    printerr("Could not malloc memory for Global->t_in_solve\n");
    exit(-1);
  } else if (Global->t_in_bar == NULL) {
    printerr("Could not malloc memory for Global->t_in_bar\n");
    exit(-1);
  } else if (Global->completion == NULL) {
    printerr("Could not malloc memory for Global->completion\n");
    exit(-1);
  }

/* POSSIBLE ENHANCEMENT:  Here is where one might distribute the a[i]
   blocks across physically distributed memories as desired.

   One way to do this is as follows:

   for (i=0;i<nblocks;i++) {
     for (j=0;j<nblocks;j++) {
       proc_num = BlockOwner(i,j);
       if ((i == nblocks-1) && (j == nblocks-1)) {
         size = edge*edge;
       } else if ((i == nblocks-1) || (j == nblocks-1)) {
         size = edge*block_size;
       } else {
         size = block_size*block_size;
       }

       Place all addresses x such that 
       (&(a[i+j*nblocks][0]) <= x < &(a[i+j*nblocks][size-1])) 
       on node proc_num
     }
   }
*/

  {
pthread_mutex_init(&((Global->start).bar_mutex), NULL);
pthread_cond_init(&((Global->start).bar_cond), NULL);
(Global->start).bar_teller=0;
};
  {pthread_mutex_init(&(Global->idlock),NULL);};
  Global->id = 0;
  for (i=1; i<P; i++) {
    {
pthread_mutex_lock(&__intern__);
assert(__threads__<__MAX_THREADS__);
pthread_create(&(__tid__[__threads__++]), NULL, (void*(*)(void *))(SlaveStart), NULL);
pthread_mutex_unlock(&__intern__);
};
  }

  InitA(rhs);
  if (doprint) {
    printf("Matrix before decomposition:\n");
    PrintA();
  }

  SlaveStart(MyNum);

  {int aantal; for(aantal=P-1;aantal>0;aantal--) pthread_join(__tid__[--__threads__], NULL);}

  if (doprint) {
    printf("\nMatrix after decomposition:\n");
    PrintA();
  }

  if (dostats) {
    maxt = avgt = mint = Global->completion[0];
    for (i=1; i<P; i++) {
      if (Global->completion[i] > maxt) {
        maxt = Global->completion[i];
      }
      if (Global->completion[i] < mint) {
        mint = Global->completion[i];
      }
      avgt += Global->completion[i];
    }
    avgt = avgt / P;
  
    min_fac = max_fac = avg_fac = Global->t_in_fac[0];
    min_solve = max_solve = avg_solve = Global->t_in_solve[0];
    min_mod = max_mod = avg_mod = Global->t_in_mod[0];
    min_bar = max_bar = avg_bar = Global->t_in_bar[0];
  
    for (i=1; i<P; i++) {
      if (Global->t_in_fac[i] > max_fac) {
        max_fac = Global->t_in_fac[i];
      }
      if (Global->t_in_fac[i] < min_fac) {
        min_fac = Global->t_in_fac[i];
      }
      if (Global->t_in_solve[i] > max_solve) {
        max_solve = Global->t_in_solve[i];
      }
      if (Global->t_in_solve[i] < min_solve) {
        min_solve = Global->t_in_solve[i];
      }
      if (Global->t_in_mod[i] > max_mod) {
        max_mod = Global->t_in_mod[i];
      }
      if (Global->t_in_mod[i] < min_mod) {
        min_mod = Global->t_in_mod[i];
      }
      if (Global->t_in_bar[i] > max_bar) {
        max_bar = Global->t_in_bar[i];
      }
      if (Global->t_in_bar[i] < min_bar) {
        min_bar = Global->t_in_bar[i];
      }
      avg_fac += Global->t_in_fac[i];
      avg_solve += Global->t_in_solve[i];
      avg_mod += Global->t_in_mod[i];
      avg_bar += Global->t_in_bar[i];
    }
    avg_fac = avg_fac/P;
    avg_solve = avg_solve/P;
    avg_mod = avg_mod/P;
    avg_bar = avg_bar/P;
  }
  printf("                            PROCESS STATISTICS\n");
  printf("              Total      Diagonal     Perimeter      Interior       Barrier\n");
  printf(" Proc         Time         Time         Time           Time          Time\n");
  printf("    0    %10.0f    %10.0f    %10.0f    %10.0f    %10.0f\n",
          Global->completion[0],Global->t_in_fac[0],
          Global->t_in_solve[0],Global->t_in_mod[0],
          Global->t_in_bar[0]);
  if (dostats) {
    for (i=1; i<P; i++) {
      printf("  %3d    %10.0f    %10.0f    %10.0f    %10.0f    %10.0f\n",
              i,Global->completion[i],Global->t_in_fac[i],
	      Global->t_in_solve[i],Global->t_in_mod[i],
	      Global->t_in_bar[i]);
    }
    printf("  Avg    %10.0f    %10.0f    %10.0f    %10.0f    %10.0f\n",
           avgt,avg_fac,avg_solve,avg_mod,avg_bar);
    printf("  Min    %10.0f    %10.0f    %10.0f    %10.0f    %10.0f\n",
           mint,min_fac,min_solve,min_mod,min_bar);
    printf("  Max    %10.0f    %10.0f    %10.0f    %10.0f    %10.0f\n",
           maxt,max_fac,max_solve,max_mod,max_bar);
  }
  printf("\n");
  Global->starttime = start;
  printf("                            TIMING INFORMATION\n");
  printf("Start time                        : %16d\n",
          Global->starttime);
  printf("Initialization finish time        : %16d\n",
          Global->rs);
  printf("Overall finish time               : %16d\n",
          Global->rf);
  printf("Total time with initialization    : %16d\n",
          Global->rf-Global->starttime);
  printf("Total time without initialization : %16d\n",
          Global->rf-Global->rs);
  printf("\n");

  if (test_result) {
    printf("                             TESTING RESULTS\n");
    CheckResult(n, a, rhs);
  }

  {exit(0);};
}


void SlaveStart()

{
  int i; 
  int j; 
  int cluster; 
  int max_block;
  int MyNum;

  {pthread_mutex_lock(&(Global->idlock));}
    MyNum = Global->id;
    Global->id ++;
  {pthread_mutex_unlock(&(Global->idlock));}

/* POSSIBLE ENHANCEMENT:  Here is where one might pin processes to
   processors to avoid migration */

  OneSolve(n, block_size, a, MyNum, dostats);
}


void OneSolve(n, block_size, a, MyNum, dostats)

double **a;
int n;
int block_size;
int MyNum;
int dostats;

{
  unsigned int i; 
  unsigned int myrs; 
  unsigned int myrf; 
  unsigned int mydone;
  struct LocalCopies *lc;

  lc = (struct LocalCopies *) malloc(sizeof(struct LocalCopies));
  if (lc == NULL) {
    fprintf(stderr,"Proc %d could not malloc memory for lc\n",MyNum);
    exit(-1);
  }
  lc->t_in_fac = 0.0;
  lc->t_in_solve = 0.0;
  lc->t_in_mod = 0.0;
  lc->t_in_bar = 0.0;

  /* barrier to ensure all initialization is done */
  {
 pthread_mutex_lock(&((Global->start).bar_mutex));
 (Global->start).bar_teller++;
 if ((Global->start).bar_teller == (P)) {
     (Global->start).bar_teller = 0;
     pthread_cond_broadcast(&((Global->start).bar_cond));
 } else
     pthread_cond_wait(&((Global->start).bar_cond), &((Global->start).bar_mutex));
 pthread_mutex_unlock(&((Global->start).bar_mutex));
};

  /* to remove cold-start misses, all processors touch their own data */
  TouchA(block_size, MyNum);

  {
 pthread_mutex_lock(&((Global->start).bar_mutex));
 (Global->start).bar_teller++;
 if ((Global->start).bar_teller == (P)) {
     (Global->start).bar_teller = 0;
     pthread_cond_broadcast(&((Global->start).bar_cond));
 } else
     pthread_cond_wait(&((Global->start).bar_cond), &((Global->start).bar_mutex));
 pthread_mutex_unlock(&((Global->start).bar_mutex));
};

/* POSSIBLE ENHANCEMENT:  Here is where one might reset the
   statistics that one is measuring about the parallel execution */

  if ((MyNum == 0) || (dostats)) {
    {long time(); (myrs) = time(0);};
  }

  lu(n, block_size, MyNum, lc, dostats);

  if ((MyNum == 0) || (dostats)) {
    {long time(); (mydone) = time(0);};
  }

  {
 pthread_mutex_lock(&((Global->start).bar_mutex));
 (Global->start).bar_teller++;
 if ((Global->start).bar_teller == (P)) {
     (Global->start).bar_teller = 0;
     pthread_cond_broadcast(&((Global->start).bar_cond));
 } else
     pthread_cond_wait(&((Global->start).bar_cond), &((Global->start).bar_mutex));
 pthread_mutex_unlock(&((Global->start).bar_mutex));
};

  if ((MyNum == 0) || (dostats)) {
    Global->t_in_fac[MyNum] = lc->t_in_fac;
    Global->t_in_solve[MyNum] = lc->t_in_solve;
    Global->t_in_mod[MyNum] = lc->t_in_mod;
    Global->t_in_bar[MyNum] = lc->t_in_bar;
    Global->completion[MyNum] = mydone-myrs;
  }
  if (MyNum == 0) {
    {long time(); (myrf) = time(0);};
    Global->rs = myrs;
    Global->done = mydone;
    Global->rf = myrf;
  }
}


void lu0(a, n, stride)

double *a;
int n; 
int stride;

{
  int j; 
  int k; 
  int length;
  double alpha;

  for (k=0; k<n; k++) {
    /* modify subsequent columns */
    for (j=k+1; j<n; j++) {
      a[k+j*stride] /= a[k+k*stride];
      alpha = -a[k+j*stride];
      length = n-k-1;
      daxpy(&a[k+1+j*stride], &a[k+1+k*stride], n-k-1, alpha);
    }
  }
}


void bdiv(a, diag, stride_a, stride_diag, dimi, dimk)

double *a; 
double *diag;
int stride_a; 
int stride_diag; 
int dimi; 
int dimk;

{
  int j; 
  int k;
  double alpha;

  for (k=0; k<dimk; k++) {
    for (j=k+1; j<dimk; j++) {
      alpha = -diag[k+j*stride_diag];
      daxpy(&a[j*stride_a], &a[k*stride_a], dimi, alpha);
    }
  }
}


void bmodd(a, c, dimi, dimj, stride_a, stride_c)

double *a; 
double *c;
int dimi; 
int dimj; 
int stride_a; 
int stride_c;

{
  int i; 
  int j; 
  int k; 
  int length;
  double alpha;

  for (k=0; k<dimi; k++) {
    for (j=0; j<dimj; j++) {
      c[k+j*stride_c] /= a[k+k*stride_a];
      alpha = -c[k+j*stride_c];
      length = dimi - k - 1;
      daxpy(&c[k+1+j*stride_c], &a[k+1+k*stride_a], dimi-k-1, alpha);
    }
  }
}


void bmod(a, b, c, dimi, dimj, dimk, stridea, strideb, stridec)

double *a; 
double *b; 
double *c;
int dimi; 
int dimj; 
int dimk; 
int stridea;
int strideb;
int stridec;

{
  int i; 
  int j; 
  int k;
  double alpha;

  for (k=0; k<dimk; k++) {
    for (j=0; j<dimj; j++) {
      alpha = -b[k+j*strideb]; 
      daxpy(&c[j*stridec], &a[k*stridea], dimi, alpha);
    }
  }
}


void daxpy(a, b, n, alpha)

double *a; 
double *b; 
double alpha;
int n;

{
  int i;

  for (i=0; i<n; i++) {
    a[i] += alpha*b[i];
  }
}


int BlockOwner(I, J)

int I; 
int J;

{
  return((J%num_cols) + (I%num_rows)*num_cols); 
}


void lu(n, bs, MyNum, lc, dostats)

int n;
int bs;
int MyNum;
struct LocalCopies *lc;
int dostats;

{
  int i, il, j, jl, k, kl;
  int I, J, K;
  double *A, *B, *C, *D;
  int dimI, dimJ, dimK;
  int strI, strJ, strK;
  unsigned int t1, t2, t3, t4, t11, t22;
  int diagowner;
  int colowner;

  for (k=0, K=0; k<n; k+=bs, K++) {
    kl = k + bs; 
    if (kl > n) {
      kl = n;
      strK = kl - k;
    } else {
      strK = bs;
    }

    if ((MyNum == 0) || (dostats)) {
      {long time(); (t1) = time(0);};
    }

    /* factor diagonal block */
    diagowner = BlockOwner(K, K);
    if (diagowner == MyNum) {
      A = a[K+K*nblocks]; 
      lu0(A, strK, strK);
    }

    if ((MyNum == 0) || (dostats)) {
      {long time(); (t11) = time(0);};
    }

    {
 pthread_mutex_lock(&((Global->start).bar_mutex));
 (Global->start).bar_teller++;
 if ((Global->start).bar_teller == (P)) {
     (Global->start).bar_teller = 0;
     pthread_cond_broadcast(&((Global->start).bar_cond));
 } else
     pthread_cond_wait(&((Global->start).bar_cond), &((Global->start).bar_mutex));
 pthread_mutex_unlock(&((Global->start).bar_mutex));
};

    if ((MyNum == 0) || (dostats)) {
      {long time(); (t2) = time(0);};
    }

    /* divide column k by diagonal block */
    D = a[K+K*nblocks];
    for (i=kl, I=K+1; i<n; i+=bs, I++) {
      if (BlockOwner(I, K) == MyNum) {  /* parcel out blocks */
	il = i + bs; 
	if (il > n) {
	  il = n;
          strI = il - i;
        } else {
          strI = bs;
        }
	A = a[I+K*nblocks]; 
	bdiv(A, D, strI, strK, strI, strK);  
      }
    }
    /* modify row k by diagonal block */
    for (j=kl, J=K+1; j<n; j+=bs, J++) {
      if (BlockOwner(K, J) == MyNum) {  /* parcel out blocks */
	jl = j+bs; 
	if (jl > n) {
	  jl = n;
          strJ = jl - j;
        } else {
          strJ = bs;
        }
        A = a[K+J*nblocks];
	bmodd(D, A, strK, strJ, strK, strK);
      }
    }

    if ((MyNum == 0) || (dostats)) {
      {long time(); (t22) = time(0);};
    }   

    {
 pthread_mutex_lock(&((Global->start).bar_mutex));
 (Global->start).bar_teller++;
 if ((Global->start).bar_teller == (P)) {
     (Global->start).bar_teller = 0;
     pthread_cond_broadcast(&((Global->start).bar_cond));
 } else
     pthread_cond_wait(&((Global->start).bar_cond), &((Global->start).bar_mutex));
 pthread_mutex_unlock(&((Global->start).bar_mutex));
};

    if ((MyNum == 0) || (dostats)) {
      {long time(); (t3) = time(0);};
    }

    /* modify subsequent block columns */
    for (i=kl, I=K+1; i<n; i+=bs, I++) {
      il = i+bs; 
      if (il > n) {
	il = n;
        strI = il - i;
      } else {
        strI = bs;
      }
      colowner = BlockOwner(I,K);
      A = a[I+K*nblocks]; 
      for (j=kl, J=K+1; j<n; j+=bs, J++) {
	jl = j + bs; 
	if (jl > n) {
	  jl = n;
          strJ= jl - j;
        } else {
          strJ = bs;
        }
	if (BlockOwner(I, J) == MyNum) {  /* parcel out blocks */
	  B = a[K+J*nblocks]; 
	  C = a[I+J*nblocks];
	  bmod(A, B, C, strI, strJ, strK, strI, strK, strI);
	}
      }
    }

    if ((MyNum == 0) || (dostats)) {
      {long time(); (t4) = time(0);};
      lc->t_in_fac += (t11-t1);
      lc->t_in_solve += (t22-t2);
      lc->t_in_mod += (t4-t3);
      lc->t_in_bar += (t2-t11) + (t3-t22);
    }
  }
}


void InitA(rhs)

double *rhs;

{
  int i, j;
  int ii, jj;
  int edge;
  int ibs;
  int jbs, skip;

  srand48((long) 1);
  edge = n%block_size;
  for (j=0; j<n; j++) {
    for (i=0; i<n; i++) {
      if ((n - i) <= edge) {
	ibs = edge;
	ibs = n-edge;
	skip = edge;
      } else {
	ibs = block_size;
	skip = block_size;
      }
      if ((n - j) <= edge) {
	jbs = edge;
	jbs = n-edge;
      } else {
	jbs = block_size;
      }
      ii = (i/block_size) + (j/block_size)*nblocks;
      jj = (i%ibs)+(j%jbs)*skip;
      a[ii][jj] = ((double) lrand48())/MAXRAND;
      if (i == j) {
	a[ii][jj] *= 10;
      }
    }
  }

  for (j=0; j<n; j++) {
    rhs[j] = 0.0;
  }
  for (j=0; j<n; j++) {
    for (i=0; i<n; i++) {
      if ((n - i) <= edge) {
	ibs = edge;
	ibs = n-edge;
	skip = edge;
      } else {
	ibs = block_size;
	skip = block_size;
      }
      if ((n - j) <= edge) {
	jbs = edge;
	jbs = n-edge;
      } else {
	jbs = block_size;
      }
      ii = (i/block_size) + (j/block_size)*nblocks;
      jj = (i%ibs)+(j%jbs)*skip;
      rhs[i] += a[ii][jj];
    }
  }
}


double TouchA(bs, MyNum)

int bs; 
int MyNum;

{
  int i, j, I, J;
  double tot = 0.0;
  int ibs;
  int jbs;

  /* touch my portion of A[] */

  for (J=0; J<nblocks; J++) {
    for (I=0; I<nblocks; I++) {
      if (BlockOwner(I, J) == MyNum) {
	if (J == nblocks-1) {
	  jbs = n%bs;
	  if (jbs == 0) {
	    jbs = bs;
          }
	} else {
	  jbs = bs;
	}
	if (I == nblocks-1) {
	  ibs = n%bs;
	  if (ibs == 0) {
	    ibs = bs;
          }
	} else {
	  ibs = bs;
	}
	for (j=0; j<jbs; j++) {
	  for (i=0; i<ibs; i++) {
	    tot += a[I+J*nblocks][i+j*ibs];
          }
	}
      }
    }
  } 
  return(tot);
}


void PrintA()

{
  int i, j;
  int ii, jj;
  int edge;
  int ibs, jbs, skip;

  edge = n%block_size;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      if ((n - i) <= edge) {
	ibs = edge;
	ibs = n-edge;
        skip = edge;
      } else {
	ibs = block_size;
        skip = block_size;
      }
      if ((n - j) <= edge) {
	jbs = edge;
	jbs = n-edge;
      } else {
	jbs = block_size;
      }
      ii = (i/block_size) + (j/block_size)*nblocks;
      jj = (i%ibs)+(j%jbs)*skip;
      printf("%8.1f ", a[ii][jj]);   
    }
    printf("\n");
  }
  fflush(stdout);
}


void CheckResult(n, a, rhs)

int n;
double **a; 
double *rhs;

{
  int i, j, bogus = 0;
  double *y, diff, max_diff;
  int ii, jj;
  int edge;
  int ibs, jbs, skip;

  edge = n%block_size;
  y = (double *) malloc(n*sizeof(double));  
  if (y == NULL) {
    printerr("Could not malloc memory for y\n");
    exit(-1);
  }
  for (j=0; j<n; j++) {
    y[j] = rhs[j];
  }
  for (j=0; j<n; j++) {
    if ((n - j) <= edge) {
      jbs = edge;
      jbs = n-edge;
      skip = edge;
    } else {
      jbs = block_size;
      skip = block_size;
    }
    ii = (j/block_size) + (j/block_size)*nblocks;
    jj = (j%jbs)+(j%jbs)*skip;

    y[j] = y[j]/a[ii][jj];
    for (i=j+1; i<n; i++) {
      if ((n - i) <= edge) {
        ibs = edge;
        ibs = n-edge;
        skip = edge;
      } else {
        ibs = block_size;
        skip = block_size;
      }
      ii = (i/block_size) + (j/block_size)*nblocks;
      jj = (i%ibs)+(j%jbs)*skip;

      y[i] -= a[ii][jj]*y[j];
    }
  }

  for (j=n-1; j>=0; j--) {
    for (i=0; i<j; i++) {
      if ((n - i) <= edge) {
	ibs = edge;
        ibs = n-edge;
        skip = edge;
      } else {
	ibs = block_size;
        skip = block_size;
      }
      if ((n - j) <= edge) {
	jbs = edge;
        jbs = n-edge;
      } else {
	jbs = block_size;
      }
      ii = (i/block_size) + (j/block_size)*nblocks;
      jj = (i%ibs)+(j%jbs)*skip;
      y[i] -= a[ii][jj]*y[j];
    }
  }

  max_diff = 0.0;
  for (j=0; j<n; j++) {
    diff = y[j] - 1.0;
    if (fabs(diff) > 0.00001) {
      bogus = 1;
      max_diff = diff;
    }
  }
  if (bogus) {
    printf("TEST FAILED: (%.5f diff)\n", max_diff);
  } else {
    printf("TEST PASSED\n");
  }
  free(y);
}


void printerr(s)

char *s;

{
  fprintf(stderr,"ERROR: %s\n",s);
}

/* Generated from ../Changed/lu.C */
