/* 
   The hpsort_index* routines are based on gsl_heapsort_index routine
   of GSL software. Here is the GSL license:
*/

/*
 * Implement Heap sort -- direct and indirect sorting
 * Based on descriptions in Sedgewick "Algorithms in C"
 * Copyright (C) 1999  Thomas Walter
 *
 * 18 February 2000: Modified for GSL by Brian Gough
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * This source is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 */

#include <stdlib.h>


/* hpsort.c */
int hpsort_index(size_t count, const double *data, int *p);
int hpsort_index1(size_t count, const double *data, int *p);
int hpsort_index_f(size_t count, const float *data, int *p);
int hpsort_index1_f(size_t count, const float *data, int *p);

#define DOWNHEAP(p, data, N, k) \
do { \
  const size_t _pki = p[k]; \
  \
  while (k <= N / 2) \
    { \
      size_t _j = 2 * k; \
      \
      if ( _j < N && data[ p[_j] ] < data[ p[_j + 1] ] ) _j++; \
      if (data[ _pki ] >= data[ p[_j] ]) break; \
      \
      p[k] = p[_j]; \
      k    = _j; \
    } \
  \
  p[k] = _pki; \
} while(0)

/*
  
DOUBLE VERSION of ROUTINES

*/


int
hpsort_index (size_t count, const double *data, int *p)
{
  /* 
     Sort the array in ascending order. This is a true inplace
     algorithm with N log N operations. Worst case (an already sorted
     array) is something like 20% slower 
  */
  
  size_t i, k, N, zero;
  
  if (count == 0) return 0; /* No data to sort */
		    
  for (i = 0; i < count; i++) p[i] = i; /* set permutation to identity */
  
  /* 
     We have n_data elements, last element is at 'n_data-1', first at
     '0' Set N to the last element number. 
  */
  
  N = count - 1;  
  k = N / 2;
  k++;  /* Compensate the first use of 'k--' */

  do
    {
      k--;
      DOWNHEAP (p, data, N, k);
    }
  while (k > 0);

  while (N > 0)
    {
      /* first swap the elements */
      size_t tmp = p[0];
      p[0] = p[N];
      p[N] = tmp;
      
      /* then process the heap */
      N--;      
      zero = 0;
      DOWNHEAP (p, data, N, zero);
    }

  return 0;
}

int
hpsort_index1 (size_t count, const double *data, int *p)
{

  /*
    Same as hpsort_index, but the array p and data 
    are in range p[1--count], and data[1--count].
    
    BEWARE: this requires the the size of p and data to be at least:
    p[count+1] && data[count+1]
  */
  
  size_t i, k, N, one;
  
  if (count < 2) {
    p[1] = 1;
    return 0; /* No data to sort */
  }
		    
  for (i = 0; i <= count; i++) p[i] = i; /* set permutation to identity */
  
  /* 
     We have n_data elements, last element is at 'n_data-1', first at
     '0' Set N to the last element number. 
  */
  
  N = count;  
  k = N / 2 + 1;
  k++;  /* Compensate the first use of 'k--' */

  do
    {
      k--;
      DOWNHEAP (p, data, N, k);
    }
  while (k > 1);

  while (N > 1)
    {
      /* first swap the elements */
      size_t tmp = p[1];
      p[1] = p[N];
      p[N] = tmp;
      
      /* then process the heap */
      N--;     
      one = 1;
      DOWNHEAP (p, data, N, one);
    }
  return 0;
}



/*
  
FLOAT VERSION of ROUTINES

*/

int
hpsort_index_f (size_t count, const float *data, int *p)
{
  /* 
     Sort the array in ascending order. This is a true inplace
     algorithm with N log N operations. Worst case (an already sorted
     array) is something like 20% slower 
  */
  
  size_t i, k, N, zero;
  
  if (count == 0) return 0; /* No data to sort */
		    
  for (i = 0; i < count; i++) p[i] = i; /* set permutation to identity */
  
  /* 
     We have n_data elements, last element is at 'n_data-1', first at
     '0' Set N to the last element number. 
  */
  
  N = count - 1;  
  k = N / 2;
  k++;  /* Compensate the first use of 'k--' */

  do
    {
      k--;
      DOWNHEAP (p, data, N, k);
    }
  while (k > 0);

  while (N > 0)
    {
      /* first swap the elements */
      size_t tmp = p[0];
      p[0] = p[N];
      p[N] = tmp;
      
      /* then process the heap */
      N--;      
      zero = 0;
      DOWNHEAP (p, data, N, zero);
    }
  return 0;
}

int
hpsort_index1_f (size_t count, const float *data, int *p)
{

  /*
    Same as hpsort_index_f, but the array p and data 
    are in range p[1--count], and data[1--count].
    
    BEWARE: this requires the the size of p and data to be at least:
    p[count+1] && data[count+1]
  */
  
  size_t i, k, N, one;
  
  if (count < 2) {
    p[1] = 1;
    return 0; /* No data to sort */
  }
		    
  for (i = 0; i <= count; i++) p[i] = i; /* set permutation to identity */
  
  /* 
     We have n_data elements, last element is at 'n_data-1', first at
     '0' Set N to the last element number. 
  */
  
  N = count;  
  k = N / 2 + 1;
  k++;  /* Compensate the first use of 'k--' */

  do
    {
      k--;
      one = 0;
      DOWNHEAP (p, data, N, k);
    }
  while (k > 1);

  while (N > 1)
    {
      /* first swap the elements */
      size_t tmp = p[1];
      p[1] = p[N];
      p[N] = tmp;
      
      /* then process the heap */
      N--;    
      one = 1;
      DOWNHEAP (p, data, N, one);
    }
  return 0;
}


/*
      int
      main()
      {
      double a[100] = {
      3000, 101, 2, 4, 5, 1, 305, 11, 58, 73, 45, 55, 77, 89, 100000, 10
      };
      int i, N = 12;
      int index[100];
      
      hpsort_index1 (N, a, index);
      for (i=1; i<N+1; i++) {
      printf("i=%d; a = %f\n", i, a[index[i]]);
      }
      }
*/
