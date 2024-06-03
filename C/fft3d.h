#ifndef FFT3D_H_
#define FFT3D_H_

#define COMPLEX_SET(dst, re, im) \
  do {				 \
    dst[0] = re;		 \
    dst[1] = im;		 \
  } while(0)

#define COMPLEX_CONJ(dst, src)   \
  do {				 \
    dst[0] =  src[0];		 \
    dst[1] = -src[1];		 \
  } while(0)

#define COMPLEX_COPY(dst, src)   \
  do {				 \
    dst[0] = src[0];		 \
    dst[1] = src[1];		 \
  } while(0)

#define COMPLEX_SCALE(comp, scale) \
  do {				 \
    comp[0] *= scale;		 \
    comp[1] *= scale;		 \
  } while(0)

#define FT_MAX(n) ( n/2 )
#define FT_MIN(n) ( (1 - n%2) - n/2 )
#define FT_INDEX(n,N) ( n<0 ? N+n : n  )

/* fft3d.c */
extern float ***general_grid_fft_interpolator_tensor3f(int n[3], int degree[3], float ***src);
extern double *fft3d_interpolator(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func);
extern void fft3d_interpolate(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func, double *result);

/* fft_util.c */
extern void fft2d_enlarge_complex_r2c(size_t n1, size_t n2, fftw_complex *in, size_t n1new, size_t n2new, fftw_complex *out);
extern void fft3d_enlarge_complex_r2c(size_t n1, size_t n2, size_t n3, fftw_complex *in, size_t n1new, size_t n2new, size_t n3new, fftw_complex *out);

#endif
