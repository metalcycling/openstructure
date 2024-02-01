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
#include "property_id.hh"

#include <boost/version.hpp>
#if BOOST_VERSION>=103800
#define BOOST_SPIRIT_USE_OLD_NAMESPACE
#include <boost/spirit/include/classic_symbols.hpp>
#else
#include <boost/spirit/symbols.hpp>
#endif


namespace ost { namespace mol {

namespace {
struct Properties : public boost::spirit::symbols<Prop> {
  Properties() {
    add
      ("rname", Prop(Prop::RNAME, Prop::STRING, Prop::RESIDUE))
      ("aname", Prop(Prop::ANAME, Prop::STRING, Prop::ATOM))  
      ("rnum", Prop(Prop::RNUM, Prop::INT, Prop::RESIDUE))
      ("rtype", Prop(Prop::RTYPE, Prop::STRING, Prop::RESIDUE))        
      ("ele", Prop(Prop::ELE, Prop::STRING, Prop::ATOM))
      ("anum", Prop(Prop::ANUM, Prop::INT, Prop::ATOM))
      ("aindex", Prop(Prop::AINDEX, Prop::INT, Prop::ATOM))
      ("x",Prop(Prop::AX, Prop::FLOAT, Prop::ATOM))
      ("y",Prop(Prop::AY, Prop::FLOAT, Prop::ATOM))
      ("z",Prop(Prop::AZ, Prop::FLOAT, Prop::ATOM))
      ("occ",Prop(Prop::OCC, Prop::FLOAT, Prop::ATOM)) 
      ("cname",Prop(Prop::CNAME, Prop::STRING, Prop::CHAIN))
      ("chain",Prop(Prop::CNAME, Prop::STRING, Prop::CHAIN))
      ("ishetatm",Prop(Prop::ISHETATM, Prop::INT, Prop::ATOM))
      ("abfac", Prop(Prop::ABFAC, Prop::FLOAT, Prop::ATOM))
      ("rbfac", Prop(Prop::RBFAC, Prop::FLOAT, Prop::RESIDUE))
      ("peptide", Prop(Prop::PEPTIDE, Prop::INT, Prop::RESIDUE))
      ("nucleotide", Prop(Prop::NUCLEOTIDE, Prop::INT, Prop::RESIDUE))
      ("rindex", Prop(Prop::RINDEX, Prop::INT, Prop::RESIDUE))
      ("protein", Prop(Prop::PROTEIN, Prop::INT, Prop::RESIDUE))
      ("ligand", Prop(Prop::LIGAND, Prop::INT, Prop::RESIDUE))
      ("water", Prop(Prop::WATER, Prop::INT, Prop::RESIDUE))
      ("acharge", Prop(Prop::ACHARGE, Prop::FLOAT, Prop::ATOM));
  }
} properties;
}

Prop::ID PropertyIDFromString(const String& prop)
{
  Prop* property=find(properties, prop.c_str());
  return NULL==property ? Prop::UNDEF : property->id;
}

Prop PropertyFromString(const String& prop)
{
  if (prop[0]=='g') {
    switch (prop[1]) {
      case 'a':
        return Prop(Prop::CUSTOM, Prop::FLOAT, Prop::ATOM);
      case 'r':
        return Prop(Prop::CUSTOM, Prop::FLOAT, Prop::RESIDUE);
      case 'c':
        return Prop(Prop::CUSTOM, Prop::FLOAT, Prop::CHAIN);
    }
  }
  Prop* property=find(properties, prop.c_str());
  return NULL==property ? Prop() : *property;
}

String Prop::GetName() const 
{
  switch (id) {
    case RNAME:
      return "rname";
    case ANAME:
      return "aname";
    case ELE:
      return "ele";
    case CNAME:
      return "cname";
    case AX:
      return "x";
    case AY:
      return "y";
    case AZ:
      return "z";
    case OCC:
      return "occ";
    case RNUM:
      return "rnum";
    case ANUM:
      return "anum";
    case ISHETATM:
      return "ishetatm";
    case RBFAC:
      return "rbfac";
    case ABFAC:
      return "abfac";
    case PEPTIDE:
      return "peptide";
    case NUCLEOTIDE:
      return "nucleotide";
    case ACHARGE:
      return "acharge";
    case RINDEX:
      return "rindex";
    case PROTEIN:
      return "protein";
    case LIGAND:
      return "ligand";
    case WATER:
      return "water";
    default:
      return "";
  }
}

String Prop::GetTypeName() const
{

  switch(type) {
    case Prop::STRING:
      return "string";
    case Prop::FLOAT:
      return "floating point";
    case Prop::INT:
      return "integral";
    case Prop::VEC_DIST:
      return "vector and distance";
  }  
  return "";
}

}} //ns

