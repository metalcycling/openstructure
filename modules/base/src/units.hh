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
#ifndef OST_UNITS_HH
#define OST_UNITS_HH

#include <ost/module_config.hh>

/*
  convenience struct for units

  Author: Ansgar Philippsen
*/

namespace ost {

/*!
  The Units struct facilitates the
  standardized usage of various units.
  It is mostly meant to ensure that all
  metric units are on the same scale; the
  standard unit in IPLT is Angstrom, which
  is defined here as 1.0. The other metric
  units are scaled accordingly, for example
  Meter is 1e10.
*/

#ifdef _MSC_VER
namespace Units {
const Real m=1.0e10;
const Real cm=1.0e8;
const Real mm=1.0e7;
const Real um=1.0e4;
const Real nm=1.0e1;
const Real pm=1.0e-2;
const Real A=1.0;

const Real in=2.54e8;
const Real ft=3.048e9;
const Real yd=9.144e9;

// standard angular unit is radians
const Real deg=0.0174532925199;
const Real rad=1.0;
const Real mrad=1.0e-3;

// energy units converted according to metric units
const Real J=1.0e20;
const Real eV=16.0217653;


//  converted according to metric units
const Real V=1.0e20;
const Real kV=1.0e23;

// time units
const Real s=1.0;
const Real min=60.0;
const Real h=3600.0;

//mass units
const Real kg=1.0;
const Real g=1.0e-3;
const Real mg=1.0e-6;
const Real kDa=1.6605655e-24;
const Real Da=1.6605655e-27;

//charge
const Real C=1.0;

} // Units

#else

struct DLLEXPORT_OST_BASE Units {
  static const Real m;
  static const Real cm;
  static const Real mm;
  static const Real um;
  static const Real nm;
  static const Real pm;
  static const Real A;

  static const Real in;
  static const Real ft;
  static const Real yd;

  static const Real deg;
  static const Real rad;
  static const Real mrad;

  static const Real J;
  static const Real eV;

  static const Real V;
  static const Real kV;

  static const Real s;
  static const Real min;
  static const Real h;

  static const Real kg;
  static const Real g;
  static const Real mg;
  static const Real kDa;
  static const Real Da;

  static const Real C;
};

#endif

} // ns







#endif // OST_UNITS_HH
