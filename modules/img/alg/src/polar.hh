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
#ifndef IMG_ALG_POLAR_HH
#define IMG_ALG_POLAR_HH

#include <ost/base.hh>
#include <ost/img/algorithm.hh>


/*
  author: Andreas Schenk
*/


namespace ost { namespace img { namespace alg {

class DLLEXPORT Polar: public ModOPAlgorithm
{
public:
    Polar();
    Polar(Real sampling);
  // algorithm interface
  ImageHandle Visit(const ConstImageHandle& i);
private:
    Real sampling_;
};

class DLLEXPORT EllipticPolar: public ModOPAlgorithm
{
public:
    EllipticPolar(Real ratio=1,Real angle=0,Real sampling=1);
  // algorithm interface
  ImageHandle Visit(const ConstImageHandle& i);
private:
    Real ratio_;
    Real angle_;
    Real sampling_;
};

class DLLEXPORT Cartesian: public ModOPAlgorithm
{
public:
    Cartesian();
  // algorithm interface
  ImageHandle Visit(const ConstImageHandle& i);
private:
};

}}} // ns

#endif // IPLT_ALG_POLAR_HH
