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
#include <ost/message.hh>
#include <ost/mol/mol.hh>
#include <ost/log.hh>

#include "conop.hh"
#include <ost/profile.hh>

namespace ost { namespace conop {

// Meyers singleton implementation
Conopology& Conopology::Instance()
{
  static Conopology impl;
  return impl;
}

Conopology::Conopology():
  known_elements_(),
  lib_()
{
  
  known_elements_.insert("AC");
  known_elements_.insert("AG");
  known_elements_.insert("AL");
  known_elements_.insert("AM");
  known_elements_.insert("AR");
  known_elements_.insert("AS");
  known_elements_.insert("AT");
  known_elements_.insert("AU");
  known_elements_.insert("B");
  known_elements_.insert("BA");
  known_elements_.insert("BE");
  known_elements_.insert("BH");
  known_elements_.insert("BI");
  known_elements_.insert("BK");
  known_elements_.insert("BR");
  known_elements_.insert("C");
  known_elements_.insert("CA");
  known_elements_.insert("D");
  known_elements_.insert("CD");
  known_elements_.insert("CE");
  known_elements_.insert("CF");
  known_elements_.insert("CL");
  known_elements_.insert("CM");
  known_elements_.insert("CO");
  known_elements_.insert("CR");
  known_elements_.insert("CS");
  known_elements_.insert("CU");
  known_elements_.insert("DB");
  known_elements_.insert("DS");
  known_elements_.insert("DY");
  known_elements_.insert("ER");
  known_elements_.insert("ES");
  known_elements_.insert("EU");
  known_elements_.insert("F");
  known_elements_.insert("FE");
  known_elements_.insert("FM");
  known_elements_.insert("FR");
  known_elements_.insert("GA");
  known_elements_.insert("GD");
  known_elements_.insert("GE");
  known_elements_.insert("H");
  known_elements_.insert("HE");
  known_elements_.insert("HF");
  known_elements_.insert("HG");
  known_elements_.insert("HO");
  known_elements_.insert("HS");
  known_elements_.insert("I");
  known_elements_.insert("IN");
  known_elements_.insert("IR");
  known_elements_.insert("K");
  known_elements_.insert("KR");
  known_elements_.insert("LA");
  known_elements_.insert("LI");
  known_elements_.insert("LR");
  known_elements_.insert("LU");
  known_elements_.insert("MD");
  known_elements_.insert("MG");
  known_elements_.insert("MN");
  known_elements_.insert("MO");
  known_elements_.insert("MT");
  known_elements_.insert("N");
  known_elements_.insert("NA");
  known_elements_.insert("NB");
  known_elements_.insert("ND");
  known_elements_.insert("NE");
  known_elements_.insert("NI");
  known_elements_.insert("NO");
  known_elements_.insert("NP");
  known_elements_.insert("O");
  known_elements_.insert("OS");
  known_elements_.insert("P");
  known_elements_.insert("PA");
  known_elements_.insert("PB");
  known_elements_.insert("PD");
  known_elements_.insert("PM");
  known_elements_.insert("PO");
  known_elements_.insert("PR");
  known_elements_.insert("PT");
  known_elements_.insert("PU");
  known_elements_.insert("RA");
  known_elements_.insert("RB");
  known_elements_.insert("RE");
  known_elements_.insert("RF");
  known_elements_.insert("RH");
  known_elements_.insert("RN");
  known_elements_.insert("RU");
  known_elements_.insert("S");
  known_elements_.insert("SB");
  known_elements_.insert("SC");
  known_elements_.insert("SE");
  known_elements_.insert("SG");
  known_elements_.insert("SI");
  known_elements_.insert("SM");
  known_elements_.insert("SN");
  known_elements_.insert("SR");
  known_elements_.insert("TA");
  known_elements_.insert("TB");
  known_elements_.insert("TC");
  known_elements_.insert("TE");
  known_elements_.insert("TH");
  known_elements_.insert("TI");
  known_elements_.insert("TL");
  known_elements_.insert("TM");
  known_elements_.insert("U");
  known_elements_.insert("V");
  known_elements_.insert("W");
  known_elements_.insert("XE");
  known_elements_.insert("Y");
  known_elements_.insert("YB");
  known_elements_.insert("ZN");
  known_elements_.insert("ZR");
}

bool Conopology::IsValidElement(const String& ele) const
{
  String upper_ele=ele;
  std::transform(upper_ele.begin(),upper_ele.end(),upper_ele.begin(),toupper);
  return known_elements_.find(upper_ele)!=known_elements_.end();
}

}} // ns
