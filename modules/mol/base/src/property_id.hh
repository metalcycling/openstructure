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
#ifndef OST_PROPERTY_ID_HH
#define OST_PROPERTY_ID_HH

/*
  Author: Marco Biasini
 */
#include <ost/message.hh>

#include <ost/mol/module_config.hh>

namespace ost { namespace mol {

struct DLLEXPORT_OST_MOL Prop {
public:
  /// \brief Defined properties.
  ///
  /// CUSTOM is for GenericFloatProperties. It must always be the last, because
  /// it is used as an index to determine in which order the generic properties
  /// where specified in the query string.
  /// GenericProperties are defined with gapropname, grpropname, gcpropname for
  /// the GenericProp 'propname' at the atom, residue and chain level
  /// respectively.
  typedef enum {
    RNAME, ANAME, CNAME, ELE, RNUM, ANUM, AX, AY, AZ, OCC, AINDEX, RTYPE,
    ISHETATM, RBFAC, ABFAC, PEPTIDE, NUCLEOTIDE, ACHARGE, RINDEX, PROTEIN,
    LIGAND, WATER, WITHIN, UNDEF, CUSTOM
  } ID;
  
  typedef enum {
    ATOM=0, RESIDUE=1, CHAIN=2, UNSPECIFIED=3
  } Level;  

  typedef enum {
    STRING, INT, FLOAT, VEC_DIST
  } Type;
public:
  ID    id;
  Type  type;
  Level level;
  
  Prop(ID i, Type t, Level l):
    id(i), type(t), level(l)
  { }
  
  Prop(): 
    id(UNDEF),  type(INT), level(ATOM)
  { }
  
  bool operator==(const Prop& rhs) const
  {
    return id==rhs.id &&  type==rhs.type && level==rhs.level;
  }
  
  bool operator!=(const Prop& rhs) const
  {
    return !this->operator==(rhs);
  }  
  
  String GetName() const;
  String GetTypeName() const;
};

/// \return the property id. Prop::UNDEF is returned in case the property is
///         not known.
Prop::ID DLLEXPORT_OST_MOL PropertyIDFromString(const String& prop);

/// \return the property. Prop::UNDEF is returned in case the property is
///         not known.
Prop DLLEXPORT_OST_MOL PropertyFromString(const String& prop);


struct DLLEXPORT PropertyError: public ost::Error
{
  PropertyError(ost::mol::Prop::ID prop):
    ost::Error("invalid property"),
    prop_(prop)
  {}
  virtual ~PropertyError() throw() {}
  
  ost::mol::Prop::ID Prop() const 
  {
    return prop_;
  }
private:
  ost::mol::Prop::ID prop_;
};

}} //ns


#endif
