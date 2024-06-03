#ifndef SHANKLAND_H_
#define SHANKLAND_H_

#include <meschach/matrix2.h>
#include <fftw3.h>
#include "fft3d.h"

extern const double two_pi;

#ifndef ABS
#   define ABS(x)   ( (x)>0 ? (x) : -(x)  )
#endif

#ifndef MAX
#   define MAX(x,y) ( (x)>(y) ? (x) : (y) )
#endif

#ifndef MIN
#   define MIN(x,y) ( (x)<(y) ? (x) : (y) )
#endif

/* shank_util.c */
extern double one_func(double m2, double sigma);
extern double first_order(double m2, double sigma);
extern double first_sqrt_order(double m2, double sigma);
extern double second_order(double m2, double sigma);
extern double pickett(double m2, double sigma);
extern double gauss(double m2, double sigma);
extern double exp_abs(double m2, double sigma);
extern double estimate_sigma(double power, int m, double (*roughnessFunc)(double m2, double sigma));
extern void print_double_vector(FILE *stream, const char *name, size_t m, double *mat);
extern void print_double_matrix(FILE *stream, const char *name, size_t m, size_t n, double *mat);
extern void print_MAT_matrix(FILE *stream, const char *name, MAT *mat);

#ifdef DEBUG
   extern void print_complex_matrix(FILE *stream, const char *name, size_t m, size_t n, fftw_complex *mat); 
   extern void print_complex_vector(FILE *stream, const char *name, size_t m, fftw_complex *mat);
#endif 


/* shank_Hinv.c */
extern MAT *shank_HinvObj(size_t Npoints, size_t Nwaves, double sigma, double (*roughnessFunc)(double m2, double sigma));
extern void shank_Hinv(size_t Npoints, size_t Nwaves, double sigma, double (*roughnessFunc)(double m2, double sigma), MAT *Hinv);
extern double *shank1d_lambdaObj(int n1, int n1new, double sigma, double *func, double (*roughnessFunc)(double m2, double sigma), double *f_aver);
extern double shank1d_lambda(int n1, int n1new, double sigma, double *func, double *lam, double (*roughnessFunc)(double m2, double sigma));
extern double *shank2d_lambdaObj(int n1, int n2, int n1new, int n2new, double sigma, double *func, double (*roughnessFunc)(double m2, double sigma), double *f_aver);
extern double shank2d_lambda(int n1, int n2, int n1new, int n2new, double sigma, double *func, double *lam, double (*roughnessFunc)(double m2, double sigma));
extern double *shank3d_lambdaObj(int n1, int n2, int n3, int n1new, int n2new, int n3new, double sigma, double *func, double (*roughnessFunc)(double m2, double sigma), double *f_aver);
extern double shank3d_lambda(int n1, int n2, int n3, int n1new, int n2new, int n3new, double sigma, double *func, double *lam, double (*roughnessFunc)(double m2, double sigma));

#endif
