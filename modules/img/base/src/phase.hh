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
  author: Andreas Schenk
*/


#ifndef PHASE_HH_
#define PHASE_HH_

#include <boost/operators.hpp>
#include <iostream>
#include <ost/img/module_config.hh>
#include <ost/img/data_types.hh>

namespace ost { namespace img {
/// \brief Manages phases og Complex numbers
///
/// Manages phases of Complex numbers, and allows sums and subtractions of phases.
/// Takes correctly into account phase periodicity
///
/// \sa \ref view_phase_diff.py "View Phase Difference Example"
class DLLEXPORT_OST_IMG_BASE Phase:
    private boost::equality_comparable<Phase>,
    private boost::additive<Phase>,
    private boost::multiplicative<Phase, Real>
{
public:
  Phase():angle_(0){};
  Phase(Real angle);
  Phase(Complex comp);
  operator Real() const;
  //! addable op
  const Phase operator+=(const Phase& rhs);
  //! subtractable op
  const Phase operator-=(const Phase& rhs);
  //! negateable
  Phase operator-() const;
  //! multipliable
  const Phase operator*=(Real d);
  //! dividable
  const Phase operator/=(Real d);
protected:
 Real reduce(Real in);
 Real angle_;
};

DLLEXPORT_OST_IMG_BASE std::ostream& operator<<(std::ostream&,const Phase a);

}}//ns

#endif /*PHASE_HH_*/

/// \example view_phase_diff.py
///
/// This scripts displays the phase difference (in degrees) between corresponding pixels in the
/// Fourier Transforms of the two input images, which must be of the same size. The Fourier
/// Transforms honor the origin of the reference system, which is assumed to be at the center of
/// the two images
///
/// Usage:
///
/// \code dng view_phase_diff.py <image1> <image2> \endcode
///
/// <BR>
/// <BR>


