//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------
/*
Author: Juergen Haas
*/
#ifndef IMG_ALG_FFTW_HELPER_HH
#define IMG_ALG_FFTW_HELPER_HH

#include <ost/config.hh>


#if OST_DOUBLE_PRECISION
#define OST_FFTW_fftw_complex fftw_complex
#define OST_FFTW_fftw_plan fftw_plan
#define OST_FFTW_fftw_execute fftw_execute
#define OST_FFTW_fftw_destroy_plan fftw_destroy_plan
#define OST_FFTW_fftw_complex fftw_complex
#define OST_FFTW_fftw_plan fftw_plan
#define OST_FFTW_fftw_execute fftw_execute
#define OST_FFTW_fftw_destroy_plan fftw_destroy_plan
#define OST_FFTW_fftw_complex fftw_complex
#define OST_FFTW_fftw_plan fftw_plan
#define OST_FFTW_fftw_execute fftw_execute
#define OST_FFTW_fftw_destroy_plan fftw_destroy_plan
#define OST_FFTW_fftw_complex fftw_complex
#define OST_FFTW_fftw_plan fftw_plan
#define OST_FFTW_fftw_execute fftw_execute
#define OST_FFTW_fftw_destroy_plan fftw_destroy_plan
#define OST_FFTW_fftw_plan_dft_r2c fftw_plan_dft_r2c
#define OST_FFTW_fftw_plan_dft fftw_plan_dft
#define OST_FFTW_fftw_plan_dft_c2r fftw_plan_dft_c2r
#define OST_FFTW_fftw_plan_many_dft fftw_plan_many_dft
#else
#define OST_FFTW_fftw_complex fftwf_complex
#define OST_FFTW_fftw_plan fftwf_plan
#define OST_FFTW_fftw_execute fftwf_execute
#define OST_FFTW_fftw_destroy_plan fftwf_destroy_plan
#define OST_FFTW_fftw_complex fftwf_complex
#define OST_FFTW_fftw_plan fftwf_plan
#define OST_FFTW_fftw_execute fftwf_execute
#define OST_FFTW_fftw_destroy_plan fftwf_destroy_plan
#define OST_FFTW_fftw_complex fftwf_complex
#define OST_FFTW_fftw_plan fftwf_plan
#define OST_FFTW_fftw_execute fftwf_execute
#define OST_FFTW_fftw_destroy_plan fftwf_destroy_plan
#define OST_FFTW_fftw_complex fftwf_complex
#define OST_FFTW_fftw_plan fftwf_plan
#define OST_FFTW_fftw_execute fftwf_execute
#define OST_FFTW_fftw_destroy_plan fftwf_destroy_plan
#define OST_FFTW_fftw_plan_dft_r2c fftwf_plan_dft_r2c
#define OST_FFTW_fftw_plan_dft fftwf_plan_dft
#define OST_FFTW_fftw_plan_dft_c2r fftwf_plan_dft_c2r
#define OST_FFTW_fftw_plan_many_dft fftwf_plan_many_dft
#endif

#if OST_FFT_USE_THREADS
  #if OST_DOUBLE_PRECISION
    #define OST_FFTW_fftw_init_threads fftw_init_threads
    #define OST_FFTW_fftw_cleanup fftw_cleanup_threads
    #define OST_FFTW_fftw_plan_with_nthreads fftw_plan_with_nthreads
  #else
    #define OST_FFTW_fftw_init_threads fftwf_init_threads
    #define OST_FFTW_fftw_cleanup fftwf_cleanup_threads
    #define OST_FFTW_fftw_plan_with_nthreads fftwf_plan_with_nthreads
  #endif
#else
  void fftw_noop(unsigned int i=0){}
  #define OST_FFTW_fftw_init_threads fftw_noop
  #define OST_FFTW_fftw_plan_with_nthreads fftw_noop
  #if OST_DOUBLE_PRECISION
    #define OST_FFTW_fftw_cleanup fftw_cleanup
  #else
    #define OST_FFTW_fftw_cleanup fftwf_cleanup
  #endif
#endif

#endif
