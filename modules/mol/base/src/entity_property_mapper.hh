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
#ifndef OST_ENTITY_PROPERTY_MAPPER_HH
#define OST_ENTITY_PROPERTY_MAPPER_HH

/*
  work in progress

  Author: Ansgar Philippsen, Marco Biasini
*/

#include <vector>
#include <ost/mol/mol.hh>
#include <ost/mol/module_config.hh>

namespace ost { namespace mol {


class DLLEXPORT_OST_MOL EntityPropertyMapper {
public:
  /// \brief construct new property mapper for the property of the given name
  /// 
  /// \param prop_name is the name of the property. The name may either be the
  ///     the name of a built-in property or refer to a 
  ///     \ref generic_properties "generic property".
  /// \param level hints at the hierarchy level that stores the property. If the
  ///      property is one of the built-in properties of atoms, residues or
  ///      chains, this parameter has no effect.
  EntityPropertyMapper(const String& prop_name, 
                       Prop::Level level=Prop::UNSPECIFIED);
  /// \brief get float property throwing an exception when the value does
  ///        not exist
  /// 
  /// Only float and int properties are considered. Properties from residues
  /// and chains may be used as well.
  Real Get(const AtomHandle& atom) const;
  /// \brief get float property, returning default if the property does not
  ///        exist
  Real Get(const AtomHandle& atom, Real default_value) const;
  
  /// \brief get float property throwing an exception when the value does
  ///        not exist
  /// 
  /// Only float and int properties are considered. Properties from residues
  /// and chains may be used as well.
  Real Get(const AtomView& atom) const;
  /// \brief get float property, returning default if the property does not
  ///        exist
  Real Get(const AtomView& atom, Real default_value) const;  
  

  /// \brief get float property throwing an exception when the value does
  ///        not exist
  ///
  /// Only float and int properties are considered. Properties from residues
  /// and chains may be used as well.
  Real Get(const ResidueHandle& atom) const;
  /// \brief get float property, returning default if the property does not
  ///        exist
  Real Get(const ResidueHandle& atom, Real default_value) const;

  /// \brief get float property throwing an exception when the value does
  ///        not exist
  ///
  /// Only float and int properties are considered. Properties from residues
  /// and chains may be used as well.
  Real Get(const ResidueView& atom) const;
  /// \brief get float property, returning default if the property does not
  ///        exist
  Real Get(const ResidueView& atom, Real default_value) const;







  /// \brief get float property throwing an exception when the value does
  ///        not exist
  ///
  /// Only float and int properties are considered. Properties from residues
  /// and chains may be used as well.
  Real Get(const ChainHandle& atom) const;
  /// \brief get float property, returning default if the property does not
  ///        exist
  Real Get(const ChainHandle& atom, Real default_value) const;

  /// \brief get float property throwing an exception when the value does
  ///        not exist
  ///
  /// Only float and int properties are considered. Properties from residues
  /// and chains may be used as well.
  Real Get(const ChainView& atom) const;
  /// \brief get float property, returning default if the property does not
  ///        exist
  Real Get(const ChainView& atom, Real default_value) const;






private:
  template <typename T, bool B>
  Real get_property(const T& atom, Real default_value) const;  
  Prop      prop_;
  String    prop_name_;
};

}} // ns

#endif
