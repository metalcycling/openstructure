//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  authors: Valerio Mariani, Andreas Schenk
*/


#ifndef IMG_ALG_FILTER_FOURIER_FILTERS_H
#define IMG_ALG_FILTER_FOURIER_FILTERS_H

#include <ost/img/algorithm.hh>
#include <ost/img/image_state.hh>
#include <ost/img/image_state.hh>

namespace ost { namespace img { namespace alg { 

//! Low Pass Filter 
/*! limit = cutoff in voxel sampling units */
class DLLEXPORT_IMG_ALG LowPassFilter : public ConstModIPAlgorithm
{
  public:

    LowPassFilter(Real limit=1.0):
      ConstModIPAlgorithm("LowPassFilter"),
      limit_(1.0/limit)
    {if(limit_==0.0)throw("Invalid limit");}
	         
    virtual void Visit(ImageHandle& ih) const ;

    void SetLimit(Real limit)       {if(limit_==0.0)throw("Invalid limit");limit_=1.0/limit;}
    Real GetLimit()                 {return 1.0/limit_;}

  private:

    Real limit_;
};

//! High Pass Filter 
/*! limit = cutoff in voxel sampling units */
class DLLEXPORT_IMG_ALG HighPassFilter : public ConstModIPAlgorithm
{
  public:

    HighPassFilter(Real limit=1.0):
      ConstModIPAlgorithm("HighPassFilter"),
      limit_(1.0/limit)
    {if(limit_==0.0)throw("Invalid limit");}
	         
    virtual void Visit(ImageHandle& ih) const ;

    void SetLimit(Real limit)       {if(limit_==0.0)throw("Invalid limit");limit_=1.0/limit;}
    Real GetLimit()                 {return 1.0/limit_;}

  private:

    Real limit_;
};

//! Gaussian Low Pass Filter 
/*! limit = cutoff in voxel sampling units */
/// \sa \ref spoke_pattern.py "Spoke Pattern Example"
class DLLEXPORT_IMG_ALG GaussianLowPassFilter : public ConstModIPAlgorithm
{
  public:

    GaussianLowPassFilter(Real limit=1.0):
      ConstModIPAlgorithm("GaussianLowPassFilter"),
      limit_(1.0/limit)
    {if(limit_==0.0)throw("Invalid limit");}
	         
    virtual void Visit(ImageHandle& ih) const ;

    void SetLimit(Real limit)       {if(limit_==0.0)throw("Invalid limit");limit_=1.0/limit;}
    Real GetLimit()                 {return 1.0/limit_;}

  private:

    Real limit_;
};

//! Gaussian High Pass Filter 
/*! limit = cutoff in voxel sampling units */
class DLLEXPORT_IMG_ALG GaussianHighPassFilter : public ConstModIPAlgorithm
{
  public:

    GaussianHighPassFilter(Real limit=1.0):
      ConstModIPAlgorithm("GaussianHighPassFilter"),
      limit_(1.0/limit)
    {if(limit_==0.0)throw("Invalid limit");}
	         
    virtual void Visit(ImageHandle& ih) const ;

    void SetLimit(Real limit)       {if(limit_==0.0)throw("Invalid limit");limit_=1.0/limit;}
    Real GetLimit()                 {return 1.0/limit_;}

  private:

    Real limit_;
};


//! Fermi Low Pass Filter 
/*! limit = cutoff in voxel sampling units 
    T = temperature parameter voxel sampling units        */
class DLLEXPORT_IMG_ALG FermiLowPassFilter : public ConstModIPAlgorithm
{
  public:

    FermiLowPassFilter(Real limit=1.0, Real t=1.0):
      ConstModIPAlgorithm("FermiLowPassFilter"),
      limit_(1.0/limit),
      t_(1.0/t)
    {if(limit_==0.0)throw("Invalid limit");if(t_==0.0)throw("Invalid t");}
	         
    virtual void Visit(ImageHandle& ih) const ;

    void SetLimit(Real limit)       {if(limit_==0.0)throw("Invalid limit");limit_=1.0/limit;}
    Real GetLimit()                 {return 1.0/limit_;}

    void SetT(Real t)               {if(t_==0.0)throw("Invalid t");t_=1.0/t;}
    Real GetT()                     {return 1.0/t_;}

  private:

    Real limit_;
    Real t_;
};

//! Fermi High Pass Filter 
/*! limit = cutoff in voxel sampling units 
    T = temperature parameter voxel sampling units        */
class DLLEXPORT_IMG_ALG FermiHighPassFilter : public ConstModIPAlgorithm
{
  public:

    FermiHighPassFilter(Real limit=1.0, Real t=1.0):
      ConstModIPAlgorithm("FermiLowPassFilter"),
      limit_(1.0/limit),
      t_(-1.0/t)
    {if(limit_==0.0)throw("Invalid limit");if(t_==0.0)throw("Invalid t");}
	         
    virtual void Visit(ImageHandle& ih) const ;

    void SetLimit(Real limit)       {if(limit_==0.0)throw("Invalid limit");limit_=1.0/limit;}
    Real GetLimit()                 {return 1.0/limit_;}

    void SetT(Real t)               {if(t_==0.0)throw("Invalid t");t_=-1.0/t;}
    Real GetT()                     {return -1.0/t_;}

  private:

    Real limit_;
    Real t_;
};

//! Butterworth Low Pass Filter 
/*! limit = passband in voxel sampling units 
    stop = stopband in voxel sampling units        */
class DLLEXPORT_IMG_ALG ButterworthLowPassFilter : public ConstModIPAlgorithm
{
  public:

    ButterworthLowPassFilter(Real limit=1.0,Real stop=1.0):
      ConstModIPAlgorithm("ButterworthLowPassFilter"),
      limit_(1.0/limit),
      stop_(1.0/stop),
      eps_(0.882), 
      a_(10.624)
    {if(limit_==0.0)throw("Invalid limit");if(stop==0.0)throw("Invalid stop");}
	         
    virtual void Visit(ImageHandle& ih) const ;

    void SetLimit(Real limit)       {if(limit_==0.0)throw("Invalid limit");limit_=1.0/limit;}
    Real GetLimit()                 {return 1.0/limit_;}

    void SetStop(Real stop) {if(stop==0.0)throw("Invalid stop");stop_=1.0/stop;}
    Real GetStop()          {return 1.0/stop_;}

    void SetEps(Real eps)   {eps_=eps;}
    Real GetEps()           {return eps_;}

    void SetA(Real a)       {a_=a;}
    Real GetA()             {return a_;}
      
  private:

    Real limit_;
    Real stop_;
    Real eps_;
    Real a_;       
};

//! Butterworth High Pass Filter 
/*! limit = passband in voxel sampling units 
    stop = stopband in voxel sampling units        */
class DLLEXPORT_IMG_ALG ButterworthHighPassFilter : public ButterworthLowPassFilter
{
public:
  ButterworthHighPassFilter(Real limit=1.0,Real stop=1.0):
    ButterworthLowPassFilter(limit,stop){}

};

}}} //ns

#endif // IPLT_ALG_FILTER_FOURIER_FILTERS_H

/// \example spoke_pattern.py
///
///This script can be used to create a "spoke pattern" image in the style of the ones featured in the following paper:
///
/// Downing K.H., Glaeser R.M., Restoration of weak phase-contrast images recorded with a high degree of
/// defocus: The "twin image" problem associated with CTF correction, Ultramicroscopy 108 (2008), 921-928
///
/// All tweakable parameters are defined in the first lines of the script, which can be edited before the script is
/// executed. The script generates the image, opens a viewer to show it and saves it in the current folder with the
/// name specified by the user.
///
/// Usage:
///
/// \code dng spoke_pattern.py <output image> \endcode
///
/// <BR>
/// <BR>
