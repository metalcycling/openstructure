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
  convenience struct for units

  Author: Ansgar Philippsen
*/

#ifndef _MSC_VER

#include <cmath>

#include "units.hh"

namespace ost {

// all metric units are given in Angstrom units!
const Real Units::m=1.0e10;
const Real Units::cm=1.0e8;
const Real Units::mm=1.0e7;
const Real Units::um=1.0e4;
const Real Units::nm=1.0e1;
const Real Units::pm=1.0e-2;
const Real Units::A=1.0;

const Real Units::in=2.54e8;
const Real Units::ft=3.048e9;
const Real Units::yd=9.144e9;

// standard angular unit is radians
const Real Units::deg=M_PI/180.0;
const Real Units::rad=1.0;
const Real Units::mrad=1.0e-3;

// energy units converted according to metric units
const Real Units::J=1.0e20;
const Real Units::eV=16.0217653;


//  converted according to metric units
const Real Units::V=1.0e20;
const Real Units::kV=1.0e23;

// time units
const Real Units::s=1.0;
const Real Units::min=60.0;
const Real Units::h=3600.0;

//mass units
const Real Units::kg=1.0;
const Real Units::g=1.0e-3;
const Real Units::mg=1.0e-6;
const Real Units::kDa=1.6605655e-24;
const Real Units::Da=1.6605655e-27;

//charge
const Real Units::C=1.0;


} // ns

#endif
