#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <fftw/fftw3.h>
#include <ipps.h>
#include "iscl_hack.h"

/*!
 * These are some hacked routines from ISTI's ISCL.  They are copyrighted
 * under the Apache 2 license.
 */

//enum isclError_enum
int fft_rfftfreqs64f_work_hack(const int n, const double dt,
                          const int lenf, double *__restrict__ freqs)
{
    double xden;
    int i, nlim;
/*
    isclReturnArrayTooSmallError("n", n, 1);
    isclReturnNullPointerError("freqs", freqs);
    isclReturnArrayTooSmallError("lenf", lenf, n/2+1);
*/
    if (dt <= 0.0)
    {
        fprintf(stderr, "%s: invalid sample spacing=%e!\n", __func__, dt);
        return -1;
/*
        isclPrintError("Invalid sample spacing=%e!\n", dt);
        if (array_zeros64f_work(lenf, freqs) != ISCL_SUCCESS)
        {   
            isclPrintError("Error zeroing out freqs\n");
        }
        return ISCL_INVALID_INPUT;
*/
    }
    // Easy case
    if (n == 1)
    {
        freqs[0] = 0.0;
        return 0; //ISCL_SUCCESS;
    }
    // Set denominator
    nlim = n/2 + 1;
    xden = 1.0/(dt*(double) n);
    #pragma omp simd
    for (i=0; i<nlim; i++)
    {
        freqs[i] = (double) i*xden;
    }
    return 0; //ISCL_SUCCESS;
}

//enum isclError_enum
int fft_irfft64z_work_hack(const int nx, const double complex *__restrict__ x,
                      const int n, double *__restrict__ y)
{ 
    fftw_plan p;
    fftw_complex *in;
    double *out, xnorm;
    int i, ntf;
    double complex zero = 0.0;
/*
    isclReturnArrayTooSmallError("n", n, 1);
    isclReturnArrayTooSmallError("nx", nx, 1);
    isclReturnNullPointerError("x", x);
    isclReturnNullPointerError("y", y);
*/
    ntf = n/2 + 1;
    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(size_t) ntf);
    out = y; //(double *)fftw_malloc(sizeof(double)*(size_t) n);
    p = fftw_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
    //  Equal size transforms
    if (nx == ntf)
    {
        for (i=0; i<nx; i++)
        {
            in[i] = x[i];
        }
    }
    // Truncate x to length of output array y
    else if (nx > ntf)
    {
        for (i=0; i<ntf; i++)
        {
            in[i] = x[i];
        }
    }
    // Pad to length of transform 
    else
    {
        for (i=0; i<nx; i++)
        {
            in[i] = x[i];
        }
        for (i=nx; i<ntf; i++)
        {
            in[i] = zero;
        }
    }
    // Transform
    fftw_execute(p);
    // Copy and normalize
    xnorm = 1.0/(double) n;
    for (i=0; i<n; i++){y[i] = y[i]*xnorm;} //cblas_dscal(n, xnorm, y, 1);
    // Clean
    fftw_destroy_plan(p);
    fftw_free(in);
    //if (!__iscl_isinit()){fftw_cleanup();}
    in = NULL;
    out = NULL;
    return 0; //ISCL_SUCCESS;
}


//enum isclError_enum
int fft_fftshift64f_work_hack(const int n, const double *__restrict__ x,
                         double *__restrict__ xshift)
{
    //double *work;
    int i1, jf, ncopy1, ncopy2;
    int ierr; //enum isclError_enum ierr;
    ierr = 0;
/*
    isclReturnArrayTooSmallError("n", n, 1);
    isclReturnNullPointerError("x", x);
    isclReturnNullPointerError("xshift", xshift);
*/
    // Handle base cases explictly
    if (n == 1)
    {
        xshift[0] = x[0];
        return 0; //ISCL_SUCCESS;
    }
    if (n == 2)
    {
        xshift[0] = x[1];
        xshift[1] = x[0];
        return 0; //ISCL_SUCCESS;
    }
    // Do the general problem
    i1 = n/2;
    if (n%2 == 1)
    {
        i1 = n/2 + 1;
        ncopy1 = n - i1; // Tail shift
        ncopy2 = i1;
        jf = i1 - 1;
    }
    else
    {
        i1 = n/2;
        ncopy1 = n/2;
        ncopy2 = n/2;
        jf = i1;
    }
    // Copy second half of x to xshift 
    //ierr = array_copy64f_work(ncopy1, &x[i1], xshift); 
    ippsCopy_64f(&x[i1], xshift, ncopy1);
    if (ierr != 0) //ISCL_SUCCESS)
    {
        fprintf(stderr, "%s: error in initial copy\n", __func__); //isclPrintError("%s", "Error in initial copy");
        return ierr;
    }
    // Copy first half of x to second half of xshift
    //ierr = array_copy64f_work(ncopy2, x, &xshift[jf]);
    ippsCopy_64f(x, &xshift[jf], ncopy2);
    if (ierr != 0) //ISCL_SUCCESS)
    {
        fprintf(stderr, "%s: error in second copy\n", __func__); //isclPrintError("%s", "Error in second copy");
        return ierr;
    }
    return ierr;
}

//============================================================================//
//enum isclError_enum
int fft_rfft64f_work_hack(const int nx, const double *__restrict__ x, const int n,
                     const int ny, double complex *__restrict__ y)
{
    fftw_complex *out;
    double *in;
    fftw_plan p;
    int i;//, ntf;
    int ierr; //enum isclError_enum ierr;
    const double zero = 0.0;
    //------------------------------------------------------------------------//
    //  
    // Size checking
    ierr = 0; //ISCL_SUCCESS;
    //ntf = n/2 + 1;
    //isclReturnArrayTooSmallError("n", n, 1);
    //isclReturnArrayTooSmallError("nx", nx, 1);
    //isclReturnArrayTooSmallError("ny", ny, ntf);
    //isclReturnNullPointerError("x", x);
    //isclReturnNullPointerError("y", y);
    // Set space and make plan
    in  = (double *)fftw_malloc(sizeof(double)*(size_t) n);
    out = y; //(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(size_t) ntf);
    p = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
    // Equal size transforms
    if (nx == n)
    {
        #pragma omp simd
        for (i=0; i<n; i++)
        {
            in[i] = x[i];
        }
    }
    // Truncate x to length of output array y
    else if (nx > n)
    {
        #pragma omp simd
        for (i=0; i<n; i++)
        {
            in[i] = x[i];
        }
    }
    // Pad x to length of output array y
    else if (nx < n)
    {
        #pragma omp simd
        for (i=0; i<nx; i++)
        {
            in[i] = x[i];
        }
        #pragma omp simd
        for (i=nx; i<n; i++)
        {
            in[i] = zero;
        }
    }
    else
    {
        fprintf(stderr, "%s: Could not classify job (nx,ny)=(%d,%d)",
                __func__, nx, ny);
        //isclPrintError("Could not classify job (nx,ny)=(%d,%d)", nx, ny);
        fftw_destroy_plan(p);
        ierr =-1; //ISCL_ALGORITHM_FAILURE;
        goto ERROR;
    }
    // Transform
    fftw_execute(p);
    // Free plan and data
    fftw_destroy_plan(p);
    fftw_free(in);
    //if (!__iscl_isinit()){fftw_cleanup();} // Assumed we'll clean this later
    // don't clean up
    in = NULL;
    out = NULL;
    return 0; //ISCL_SUCCESS;
    // In case of error set output to zero
ERROR:;
/*
    if (array_zeros64z_work(ny, y) != ISCL_SUCCESS)
    {
        isclPrintError("%s", "Error nullying out y");
    }
*/
    return ierr;
}

