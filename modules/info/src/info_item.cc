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
  high level info item
  
  Author: Ansgar Philippsen
*/

#include <iostream>

#include <boost/regex.hpp>

#include "info.hh"
#include "item_type.hh"
#include "info_impl.hh"

namespace ost { namespace info {

InfoItem::InfoItem(const InfoHandle& root, const ElePtr& impl, Type t):
  root_(root),
  impl_(impl)
{
  impl_->SetStringRepr(impl_->GetAttribute("value"),true);
  impl_->SetNewType(t);
}

InfoGroup InfoItem::GetParent() const
{
  return InfoGroup(root_,impl_->GetParent());
}

String InfoItem::GetName() const
{
  return impl_->GetName();
}

void InfoItem::SetName(const String& n)
{
  impl_->SetName(n);
}


String InfoItem::GetValue() const
{
  // assume that value attribute exists
  return impl_->GetAttribute("value");
}

void InfoItem::SetValue(const String& v)
{
  impl_->SetStringRepr(v,true);
  impl_->SetType(IT_STRING);
}

Type InfoItem::GetType() const
{
  return impl_->GetType();
}

void InfoItem::SetType(Type t)
{
  impl_->SetNewType(t);
}

void InfoItem::SetType(int t)
{
  impl_->SetNewType(t);
}

#define INFO_ITEM_TYPE_WARN(name,type) std::cerr << "Warning: InfoItem '" << name << "' does not have requested type '"<< type << "'" << std::endl;

Real InfoItem::AsFloat() const
{
  if(impl_->GetType()!=IT_FLOAT) {
    INFO_ITEM_TYPE_WARN(impl_->GetName(),"float");
  }
  return impl_->GetFloatRepr();
}

int InfoItem::AsInt() const
{
  if(impl_->GetType()!=IT_INT) {
    INFO_ITEM_TYPE_WARN(impl_->GetName(),"int");
  }
  return impl_->GetIntRepr();
}

bool InfoItem::AsBool() const
{
  if(impl_->GetType()!=IT_BOOL) {
    INFO_ITEM_TYPE_WARN(impl_->GetName(),"bool");
  }
  return impl_->GetBoolRepr();
}

geom::Vec3 InfoItem::AsVector() const
{
  if(impl_->GetType()!=IT_VECTOR) {
    INFO_ITEM_TYPE_WARN(impl_->GetName(),"vec3");
  }
  return impl_->GetVecRepr();
}

#undef INFO_ITEM_TYPE_WARN

void InfoItem::SetFloat(Real d)
{
  impl_->SetFloatRepr(d);
}

void InfoItem::SetInt(int i)
{
  impl_->SetIntRepr(i);
}

void InfoItem::SetBool(bool b)
{
  impl_->SetBoolRepr(b);
}

void InfoItem::SetVector(const geom::Vec3& v)
{
  impl_->SetVecRepr(v);
}


String InfoItem::GetAttribute(const String& name) const
{
  return impl_->GetAttribute(name);
}

void InfoItem::SetAttribute(const String& name, const String& value)
{
  impl_->SetAttribute(name,value);
}

bool InfoItem::HasAttribute(const String& name) const
{
  return impl_->HasAttribute(name);
}

void InfoItem::RemoveAttribute(const String& name)
{
  impl_->RemoveAttribute(name);
}

std::vector<String> InfoItem::GetAttributeList() const
{
  return impl_->GetAttributeList();
}

bool InfoItem::operator==(const InfoItem& rhs) const
{
  return impl_==rhs.impl_;
}

bool InfoItem::operator!=(const InfoItem& rhs) const
{
  return !this->operator==(rhs);
}

String InfoItem::GetComment() const
{
  return impl_->GetComment();
}

Real GetFloatInfoItem(const InfoGroup& ig, const InfoPath& path, Real def)
{
  if(ig.HasItem(path)) {
    return ig.GetItem(path).AsFloat();
  }
  return def;
}

int GetIntInfoItem(const InfoGroup& ig, const InfoPath& path, int def)
{
  if(ig.HasItem(path)) {
    return ig.GetItem(path).AsInt();
  }
  return def;
}

bool GetBoolInfoItem(const InfoGroup& ig, const InfoPath& path, bool def)
{
  if(ig.HasItem(path)) {
    return ig.GetItem(path).AsBool();
  }
  return def;
}

geom::Vec3 GetVectorInfoItem(const InfoGroup& ig, const InfoPath& path, const geom::Vec3& def)
{
  if(ig.HasItem(path)) {
    return ig.GetItem(path).AsVector();
  }
  return def;
}

String GetStringInfoItem(const InfoGroup& ig, const InfoPath& path, const String& def)
{
  if(ig.HasItem(path)) {
    return ig.GetItem(path).GetValue();
  }
  return def;
}

void SetFloatInfoItem(InfoGroup ig, const InfoPath& path, Real val)
{
  ig.RetrieveItem(path,false).SetFloat(val);
}

void SetIntInfoItem(InfoGroup ig, const InfoPath& path, int val)
{
  ig.RetrieveItem(path,false).SetInt(val);
}

void SetBoolInfoItem(InfoGroup ig, const InfoPath& path, bool val)
{
  ig.RetrieveItem(path,false).SetBool(val);
}

void SetVectorInfoItem(InfoGroup ig, const InfoPath& path, const geom::Vec3& val)
{
  ig.RetrieveItem(path,false).SetVector(val);
}

void SetStringInfoItem(InfoGroup ig, const InfoPath& path, const String& val)
{
  ig.RetrieveItem(path,false).SetValue(val);
}

}} // ns
