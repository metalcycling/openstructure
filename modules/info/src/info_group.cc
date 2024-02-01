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

#include "info.hh"
#include "info_impl.hh"
#include "info_visitor.hh"
#include "info_error.hh"

namespace ost { namespace info {


namespace {

Type ResolveItemType(const ElePtr& p) 
{
  Type tp=IT_STRING;
  if(p->HasAttribute("type")) {
    String tps=p->GetAttribute("type");
    if(tps=="int" || tps=="INT") {
      tp=IT_INT;
    } else if(tps=="float" || tps=="FLOAT") {
      tp=IT_FLOAT;
    } else if(tps=="bool" || tps=="BOOL") {
      tp=IT_BOOL;
    } else if(tps=="vector" || tps=="VECTOR") {
	tp=IT_VECTOR;
    } else {
      // keep String default, alternatively could throw error here...
    }
  }
  return tp;
}

}


InfoGroup::InfoGroup(const InfoHandle& root, const ElePtr& impl):
  root_(root),
  impl_(impl)
{
}

InfoGroup InfoGroup::GetParent() const
{
  return InfoGroup(root_,impl_->GetParent());
}

String InfoGroup::GetName() const
{
  return impl_->GetName();
}

void InfoGroup::SetName(const String& name)
{
  impl_->SetName(name);
}

InfoPath InfoGroup::GetPath() const
{
  return InfoPath(impl_->GetPath());
}


InfoGroupList InfoGroup::GetGroups(const InfoPath& path) const
{
  std::vector<String> plist=path.GetList();
  InfoGroupList  groups;
  if(plist.empty()) {
    groups.push_back(*this);
    return groups;
  }
  std::vector<String>::const_iterator pos=plist.begin();
  this->do_group_lookup(pos, plist.end(), groups);
  return groups;
}


InfoItemList InfoGroup::GetItems(const InfoPath& path) const
{
  InfoGroupList groups=this->GetGroups(path.Strip());
  InfoItemList items;
  String iname=path.Leaf();
  for (InfoGroupList::iterator i=groups.begin(), e=groups.end(); i!=e; ++i) {
    InfoGroup group=*i;
    detail::EleList eles=group.impl_->GetSubs(iname);
    for (detail::EleList::const_iterator j=eles.begin(), 
         e2=eles.end(); j!=e2; ++j) {
       ElePtr p=*j;
       if(!p->HasAttribute("value")) {
         throw InfoError(String("requested item '") + iname + 
                         String("' has not value attribute"));
       }
       items.push_back(InfoItem(root_, *j, ResolveItemType(p)));
    }
  }
  return items;
}


InfoGroup InfoGroup::GetGroup(const InfoPath& path, bool use_defaults) const
{
  std::vector<String> plist=path.GetList();
  if(plist.empty()) return *this;
  std::vector<String>::const_iterator pos=plist.begin();
  InfoGroup last_group(*this);
  if(do_group_lookup(pos,plist.end(),last_group)) {
    return last_group;
  } else {
    InfoPath basepath=GetPath();
    if(root_.HasDefaultGroup(basepath+path) && use_defaults) {
      return root_.GetDefaultGroup(basepath+path);
    }
  }
  std::ostringstream msg;
  msg << path << " not found in group";
  throw InfoError(msg.str());
}

InfoGroup InfoGroup::CreateGroup(const String& name)
{
  return InfoGroup(root_,impl_->CreateSub(name));
}

bool InfoGroup::do_group_lookup(std::vector<String>::const_iterator& pos,
                                const std::vector<String>::const_iterator& end,
                                InfoGroupList& subgroups) const
{
  bool found=false;
  if(pos==end) {
    subgroups.push_back(InfoGroup(*this));
    return true;
  }
  detail::EleList eles=impl_->GetSubs(*pos);
  ++pos;
  if (pos==end) {
    for (detail::EleList::iterator i=eles.begin(), e=eles.end(); i!=e; ++i) {
      subgroups.push_back(InfoGroup(root_, *i));
    }
    found=!subgroups.empty();
  } else {
    for (detail::EleList::iterator i=eles.begin(), e=eles.end(); i!=e; ++i) {
      found|=InfoGroup(root_, *i).do_group_lookup(pos, end, subgroups);
    }
  }
  return found;
} 

bool InfoGroup::do_group_lookup(std::vector<String>::const_iterator& pos,
				const std::vector<String>::const_iterator& end,
				InfoGroup& subgroup) const
{
  if(pos==end) {
    return true;
  }
  if(impl_->HasSub(*pos)) {
    InfoGroup subgrp(root_,impl_->GetSub(*pos));
    // remember last valid group
    subgroup = subgrp;
    ++pos;
    if(pos==end) {
      return true;
    } else {
      return subgrp.do_group_lookup(pos,end,subgroup);
    }
  }
  return false;
} 

bool InfoGroup::HasGroup(const InfoPath& path, bool use_defaults) const
{
  std::vector<String> plist=path.GetList();
  if(plist.empty()) return true;
  std::vector<String>::const_iterator pos=plist.begin();
  InfoGroup last_group(*this);
  if(do_group_lookup(pos,plist.end(),last_group)) {
    return true;
  } else if(use_defaults){
    return root_.HasDefaultGroup(GetPath()+path);
  }
  return false; // to make compiler happy
}

InfoGroup InfoGroup::group_create(InfoGroup group, 
                                  std::vector<String>::const_iterator& pos,
                                  const std::vector<String>::const_iterator& end)
{
  InfoGroup subgroup = group.CreateGroup((*pos));
  pos++;
  if(pos==end) return subgroup;
  return group_create(subgroup,pos,end);
}

InfoGroup InfoGroup::RetrieveGroup(const InfoPath& path, bool use_defaults)
{
  if(HasGroup(path,use_defaults)) {
    return GetGroup(path);
  } else {
    std::vector<String> plist=path.GetList();
    if(plist.empty()) return *this;
    InfoGroup last_group(*this);
    std::vector<String>::const_iterator pos=plist.begin();
    bool ret=do_group_lookup(pos,plist.end(),last_group);
    if(ret) {
      return last_group;
    } else {
      return group_create(last_group,pos,plist.end());
    }
  }
}

InfoItem InfoGroup::GetItem(const InfoPath& path, bool use_defaults) const
{
  InfoGroup grp=*this;
  InfoPath grppath=path.Strip();
  std::vector<String> plist=grppath.GetList();
  std::vector<String>::const_iterator pos=plist.begin();
  String iname=path.GetList().back();
  InfoGroup last_group(*this);
  if(do_group_lookup(pos,plist.end(),last_group)) {
    if(last_group.impl_->HasSub(iname)) {

      ElePtr p(last_group.impl_->GetSub(iname));
      if(!p->HasAttribute("value")) {
        throw InfoError(String("requested item '") + iname + String("' has not value attribute"));
      }
      return InfoItem(root_,p,ResolveItemType(p));
    }
  }
  // check default
  InfoPath basepath=GetPath();
  if(use_defaults && root_.HasDefaultItem(basepath+path)) {
    return root_.GetDefaultItem(basepath+path);
  }
  throw InfoError(String("requested item '") + iname + String("' not found"));
}

InfoItem InfoGroup::CreateItem(const String& name, const String& value)
{
  ElePtr p(impl_->CreateSub(name));
  p->SetAttribute("value",value);
  return InfoItem(root_,p);
}

InfoItem InfoGroup::CreateItem(const String& name, Real value)
{
  InfoItem item=this->CreateItem(name, String(""));
  item.SetFloat(value);
  return item;
}

InfoItem InfoGroup::CreateItem(const String& name, bool value)
{
  InfoItem item=this->CreateItem(name, String(""));
  item.SetBool(value);
  return item;  
}

InfoItem InfoGroup::CreateItem(const String& name, int value)
{
  InfoItem item=this->CreateItem(name, String(""));
  item.SetInt(value);
  return item;  
}

InfoItem InfoGroup::CreateItem(const String& name, const geom::Vec3& vector)
{
  InfoItem item=this->CreateItem(name, String(""));
  item.SetVector(vector);
  return item;
}


bool InfoGroup::HasItem(const InfoPath& path, bool use_defaults) const
{
  InfoGroup grp=*this;
  InfoPath grppath=path.Strip();
  std::vector<String> plist=grppath.GetList();
  std::vector<String>::const_iterator pos=plist.begin();
  InfoGroup last_group(*this);
  if(do_group_lookup(pos,plist.end(),last_group)) {
    String iname=path.GetList().back();
    if(last_group.impl_->HasSub(iname)) {
      ElePtr p(last_group.impl_->GetSub(iname));
      if(p->HasAttribute("value")) {
        return true;
      }
    }
  } 
  // fallback to default item
  if(use_defaults){
    return root_.HasDefaultItem(GetPath()+path);
  }
  return false;
}


InfoItem InfoGroup::RetrieveItem(const InfoPath& path, bool use_defaults)
{
  String item=path.GetList().back();

  InfoGroup grp = RetrieveGroup(path.Strip(),use_defaults);
  
  if(grp.HasItem(item,use_defaults)) {
    return grp.GetItem(item);
  } else {
    return grp.CreateItem(item,"");
  }
}

void InfoGroup::Remove(const InfoPath& path, bool remove_defaults)
{
  if(HasGroup(path,remove_defaults)){
    InfoGroup g=GetGroup(path,remove_defaults);
    InfoGroup p=g.GetParent();
    p.impl_->RemoveSub(g.GetName());
  }
}

void InfoGroup::Remove(const InfoGroup& group)
{
	InfoGroup p=group.GetParent();
	p.impl_->RemoveSub(group.impl_);
}

String InfoGroup::GetAttribute(const String& name) const
{
  return impl_->GetAttribute(name);
}

void InfoGroup::SetAttribute(const String& name, const String& value)
{
  impl_->SetAttribute(name,value);
}

bool InfoGroup::HasAttribute(const String& name) const
{
  return impl_->HasAttribute(name);
}

void InfoGroup::RemoveAttribute(const String& name)
{
  impl_->RemoveAttribute(name);
}

void InfoGroup::Apply(InfoVisitor& v, bool visit_this)
{
  bool stat = visit_this ? v.VisitGroup(*this) : true;
  if(stat) {
    // descend into subgroups
    std::vector<ElePtr> epl = impl_->GetSubList();
    for(std::vector<ElePtr>::iterator it=epl.begin();it!=epl.end();++it) {
      if((*it)->HasAttribute("value")) {
        InfoItem subitem(root_,*it,ResolveItemType(*it));
        v.VisitItem(subitem);
      } else {
        InfoGroup subgroup(root_,*it);
        subgroup.Apply(v,true);
      }
    }
  }
  v.VisitGroupFinish(*this);
}

void InfoGroup::Apply(InfoConstVisitor& v, bool visit_this) const
{
  bool stat = visit_this ? v.VisitGroup(*this) : true;
  if(stat) {
    // descend into subgroups
    std::vector<ElePtr> epl = impl_->GetSubList();
    for(std::vector<ElePtr>::iterator it=epl.begin();it!=epl.end();++it) {
      if((*it)->HasAttribute("value")) {
        InfoItem subitem(root_,*it,ResolveItemType(*it));
        v.VisitItem(subitem);
      } else {
        InfoGroup subgroup(root_,*it);
        subgroup.Apply(v,true);
      }
    }
  }
  v.VisitGroupFinish(*this);
}

std::vector<String> InfoGroup::GetAttributeList() const
{
  return impl_->GetAttributeList();
}

String InfoGroup::GetTextData() const
{
  return impl_->GetTextData();
}

void InfoGroup::SetTextData(const String& td)
{
  impl_->SetTextData(td);
}

String InfoGroup::GetComment() const
{
  return impl_->GetComment();
}

bool InfoGroup::operator==(const InfoGroup& ref) const
{
  return *impl_ == *ref.impl_;
}

void InfoGroup::Swap(InfoGroup& group)
{
  std::swap(group.impl_, impl_);
  std::swap(group.root_, root_);
}

bool InfoGroup::operator!=(const InfoGroup& ref) const
{
  return !this->operator==(ref);
}

}} // ns

