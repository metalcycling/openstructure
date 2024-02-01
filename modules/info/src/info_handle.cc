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
  high level info interface
  
  Author: Ansgar Philippsen
*/

#include <iostream>
#include <sstream>

#include "info.hh"
#include "info_impl.hh"
#include "info_error.hh"

namespace ost { namespace info {

InfoHandle CreateInfo()
{
  return InfoHandle(RootPtr(new detail::InfoImpl()));
}
InfoHandle CreateInfo(const String& text)
{
  return InfoHandle(RootPtr(new detail::InfoImpl(true,text)));
}

InfoHandle LoadInfo(const String& file)
{
  return InfoHandle(RootPtr(new detail::InfoImpl(file)));
}


// handle
InfoHandle::InfoHandle():
  impl_()
{}

InfoHandle::InfoHandle(RootPtr impl):
  impl_(impl)
{}

InfoHandle InfoHandle::Copy() const
{
  return InfoHandle(impl_->Copy());
}

void InfoHandle::Import(const String& file)
{
  impl_->Import(file);
}

void InfoHandle::Export(const String& file) const
{
  impl_->Export(file);
}

bool InfoHandle::IsValid() const
{
  return static_cast<bool>(impl_);
}


InfoGroup InfoHandle::Root() const
{
  if(!impl_) throw InfoError("Empty InfoHandle");
  return InfoGroup(*this,ElePtr(new detail::EleImpl((*impl_).Root())));
}

void InfoHandle::AddDefault(const InfoHandle& h)
{
  impl_->DefList().push_back(h.impl_);
}

bool InfoHandle::HasDefaultGroup(const InfoPath& p) const
{
  for(RootPtrList::const_iterator it=impl_->DefList().begin();it!=impl_->DefList().end();++it) {
    InfoGroup lroot(*it,ElePtr(new detail::EleImpl((*it)->Root())));
    if(lroot.HasGroup(p)) {
      return true;
    }
  }
  return false;
}

InfoGroup InfoHandle::GetDefaultGroup(const InfoPath& p) const
{
  for (RootPtrList::const_iterator it=impl_->DefList().begin();it!=impl_->DefList().end();++it) {
    InfoGroup lroot(*it,ElePtr(new detail::EleImpl((*it)->Root())));
    if(lroot.HasGroup(p)) {
      return lroot.GetGroup(p);
    }
  }
  throw InfoError("default group not found");
}

bool InfoHandle::HasDefaultItem(const InfoPath& p) const
{
  if(impl_->DefList().empty()) return false;
  for(RootPtrList::const_iterator it=impl_->DefList().begin();it!=impl_->DefList().end();++it) {
    InfoGroup lroot(*it,ElePtr(new detail::EleImpl((*it)->Root())));
    if(lroot.HasItem(p)) {
      return true;
    }
  }
  return false;
}

InfoItem InfoHandle::GetDefaultItem(const InfoPath& p) const
{
  for(RootPtrList::const_iterator it=impl_->DefList().begin();it!=impl_->DefList().end();++it) {
    InfoGroup lroot(*it,ElePtr(new detail::EleImpl((*it)->Root())));
    if(lroot.HasItem(p)) {
      return lroot.GetItem(p);
    }
  }
  throw InfoError("default item not found");
}

InfoGroup InfoHandle::GetParent()  const
{
  return Root().GetParent();
}

void InfoHandle::SetName(const String& name)
{
  Root().SetName(name);
}

String InfoHandle::GetName()  const
{
  return Root().GetName();
}

InfoGroup InfoHandle::GetGroup(const InfoPath& path, bool use_defaults)  const
{
  return Root().GetGroup(path,use_defaults);
}

InfoGroup InfoHandle::CreateGroup(const String& name)
{
  return Root().CreateGroup(name);
}

bool InfoHandle::HasGroup(const InfoPath& name, bool use_defaults)  const
{
  return Root().HasGroup(name,use_defaults);
}

InfoGroup InfoHandle::RetrieveGroup(const InfoPath& path, bool use_defaults)
{
  return Root().RetrieveGroup(path,use_defaults);
}

InfoItem InfoHandle::GetItem(const InfoPath& path, bool use_defaults)  const
{
  return Root().GetItem(path,use_defaults);
}

InfoItem InfoHandle::CreateItem(const String& name, const String& value)
{
  return Root().CreateItem(name,value);
}

bool InfoHandle::HasItem(const InfoPath& path, bool use_defaults)  const
{
  return Root().HasItem(path,use_defaults);
}

InfoItem InfoHandle::RetrieveItem(const InfoPath& path, bool use_defaults)
{
  return Root().RetrieveItem(path,use_defaults);
}

void InfoHandle::Remove(const InfoPath& path, bool remove_defaults)
{
  return Root().Remove(path,remove_defaults);
}
void InfoHandle::Remove(const InfoGroup& group)
{
  return Root().Remove(group);
}

String InfoHandle::GetAttribute(const String& name)  const
{
  return Root().GetAttribute(name);
}

void InfoHandle::SetAttribute(const String& name, const String& value)
{
  Root().SetAttribute(name,value);
}

bool InfoHandle::HasAttribute(const String& name)  const
{
  return Root().HasAttribute(name);
}

void InfoHandle::RemoveAttribute(const String& name)
{
  Root().RemoveAttribute(name);
}

std::vector<String> InfoHandle::GetAttributeList()  const
{
  return Root().GetAttributeList();
}

String InfoHandle::GetTextData()  const
{
  return Root().GetTextData();
}

void InfoHandle::SetTextData(const String& td)
{
  Root().SetTextData(td);
}

void InfoHandle::Apply(InfoVisitor& v, bool visit_this)
{
  Root().Apply(v,visit_this);
}

void InfoHandle::Apply(InfoConstVisitor& v, bool visit_this)  const
{
  Root().Apply(v,visit_this);
}


}} // ns
