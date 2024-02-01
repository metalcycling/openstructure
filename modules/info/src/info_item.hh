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
  high level info item interface
  
  Author: Ansgar Philippsen
*/

#ifndef OST_DATA_INFO_ITEM_H
#define OST_DATA_INFO_ITEM_H

#include <boost/lexical_cast.hpp>

#include <ost/geom/geom.hh>
#include <ost/base.hh>

#include "info_fw.hh"
#include "info_handle.hh"
#include "info_group.hh"
#include "info_impl_fw.hh"
#include "item_type.hh"

namespace ost { namespace info {

class DLLEXPORT_OST_INFO InfoItem {
  friend class InfoGroup;
public:


  //! retrieve parent info group
  InfoGroup GetParent() const;

  //! retrieve item name
  String GetName() const;
  //! rename item
  void SetName(const String& n);

  //! retrieve value of item
  /*!
    returns the String representation of the value
  */
  String GetValue() const;

  //! set value of item
  /*!
    this sets the value as a String
  */
  void SetValue(const String& val);

  //! retrieve type
  /*!
    type refers to the encoding of the item, which is one
    of STRING (default), INT, FLOAT, BOOL or VECTOR (see InfoItem::Type)
  */
  Type GetType() const;

  //! force type
  /*
    will do internal conversion
  */
  void SetType(Type t);

  void SetType(int t);

  //! return item value as floating point
  Real AsFloat() const;
  //! return item value as integer
  int AsInt() const;
  //! return item value as boolean
  bool AsBool() const;
  //! return item value as 3d vector
  geom::Vec3 AsVector() const;

  //! set item value as floating point
  void SetFloat(Real d);
  //! set item value as integer
  void SetInt(int i);
  //! set item value as bool
  void SetBool(bool b);
  //! set item value as 3d vector
  void SetVector(const geom::Vec3& v);

  //! return attribute of given name
  String GetAttribute(const String& name) const;
  //! set attribute of given name, will be created if it does not yet exist
  void SetAttribute(const String& name, const String& value);
  //! returns true if attribute of this names exists
  bool HasAttribute(const String& name) const;
  //! remove attribute of given name
  void RemoveAttribute(const String& name);
  //! get comment just above item
  String GetComment() const;

  //! returns a list with all attribute names
  std::vector<String> GetAttributeList() const;

  bool operator==(const InfoItem& rhs) const;
  bool operator!=(const InfoItem& rhs) const;  
private:
  InfoItem(const InfoHandle& root, const ElePtr& impl, Type t=IT_STRING);

  InfoHandle root_;
  ElePtr impl_;
  Type type_;
};

//! Convenience function to extract a value from an item
template <typename T>
T Extract(const info::InfoGroup& g, const String& n)
{
  return boost::lexical_cast<T>(g.GetItem(n).GetValue());
}

// more convenience functions
Real  DLLEXPORT_OST_INFO GetFloatInfoItem(const InfoGroup& ig, const InfoPath& path, Real def);
int DLLEXPORT_OST_INFO GetIntInfoItem(const InfoGroup& ig, const InfoPath& path, int def);
bool DLLEXPORT_OST_INFO GetBoolInfoItem(const InfoGroup& ig, const InfoPath& path, bool def);
geom::Vec3 DLLEXPORT_OST_INFO GetVectorInfoItem(const InfoGroup& ig, const InfoPath& path, const geom::Vec3& def);
String DLLEXPORT_OST_INFO GetStringInfoItem(const InfoGroup& ig, const InfoPath& path, const String& def);

void DLLEXPORT_OST_INFO SetFloatInfoItem(InfoGroup ig, const InfoPath& path, Real val);
void DLLEXPORT_OST_INFO SetIntInfoItem(InfoGroup ig, const InfoPath& path, int val);
void DLLEXPORT_OST_INFO SetBoolInfoItem(InfoGroup ig, const InfoPath& path, bool val);
void DLLEXPORT_OST_INFO SetVectorInfoItem(InfoGroup ig, const InfoPath& path, const geom::Vec3& val);
void DLLEXPORT_OST_INFO SetStringInfoItem(InfoGroup ig, const InfoPath& path, const String& val);

}} // ns

#endif
