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
  low level info implementation
  
  Author: Ansgar Philippsen
*/

#ifndef OST_DATA_INFO_IMPL_H
#define OST_DATA_INFO_IMPL_H

#include <ost/message.hh>
#include <utility>
#include <vector>
#include <iostream>

#include <QtXml/QDomElement>
#include <QtXml/QDomDocument>
#include <QtXml/QDomAttr>
#include <QtXml/QDomText>

#include <ost/base.hh>
#include <ost/geom/geom.hh>

#include <ost/ptr_observer.hh>

#include "info_impl_fw.hh"
#include "item_type.hh"

namespace ost { namespace info { namespace detail {

typedef std::vector<ElePtr> EleList;
/*
  wraps a DomElement, to be used by the high level item and group class
*/
class EleImpl {
public:
  EleImpl(const QDomElement& e);

  // tag name access
  String GetName() const;
  void SetName(const String& n);

  // get full path
  String GetPath() const;


  // attribute access
  void SetAttribute(const String& name, const String& value);
  String GetAttribute(const String& name) const;
  bool HasAttribute(const String& name) const;
  void RemoveAttribute(const String& name);

  std::vector<String> GetAttributeList() const;

  String GetTextData() const;
  void SetTextData(const String& td);

  String GetComment() const;

  /// \brief get sub elements with specified element name
  /// 
  /// does not throw. returns empty list if not elements have been found
  EleList GetSubs(const String& name) const;
  
  /// \brief get single sub elements with specified name
  /// throws if no element is found
  ElePtr GetSub(const String& name) const;  
  ElePtr CreateSub(const String& name);
  bool HasSub(const String& name) const;

  void RemoveSub(const String& name);

  void RemoveSub(const ElePtr& ele);

  QDomElement& GetElement();
  const QDomElement& GetElement() const;

  InfoImpl& GetRoot();
  InfoImpl& GetRoot() const;

  std::vector<ElePtr> GetSubList() const;

  ElePtr GetParent() const;

  bool operator==(const EleImpl& ref) const;

  Type GetType() const;
  void SetType(Type t);

  void SetNewType(int t);

  void SetStringRepr(const String& v, bool settype);
  String GetStringRepr() const;
  void SetIntRepr(int repr);
  int GetIntRepr() const;
  void SetFloatRepr(float repr);
  float GetFloatRepr() const;
  void SetBoolRepr(bool repr);
  bool GetBoolRepr() const;
  void SetVecRepr(const geom::Vec3& repr);
  geom::Vec3 GetVecRepr() const;

private:
  QDomElement ele_;

  Type type_;

  String String_repr_;
  int int_repr_;
  float float_repr_;
  bool bool_repr_;
  geom::Vec3 vec_repr_;
};


class InfoImpl {
public:
  InfoImpl();
  InfoImpl(const String& file);
  InfoImpl(bool dummy,const String& text);

  RootPtr Copy() const;

  void Import(const String& file);
  void Export(const String& file);

  EleImpl Root();

  RootPtrList& DefList() {return def_list_;}
  const RootPtrList& DefList() const {return def_list_;}

  QDomDocument& GetDocument() {return document_;}
  const QDomDocument& GetDocument() const {return document_;}

private:
  QDomDocument document_;

  RootPtrList def_list_;
};

}}} // ns

#endif
