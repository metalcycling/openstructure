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
  info high level group interface
  
  Author: Ansgar Philippsen
*/

#ifndef OST_INFO_GROUP_H
#define OST_INFO_GROUP_H

#include <vector>
#include <ost/base.hh>
#include <ost/geom/geom.hh>
#include "info_fw.hh"
#include "info_handle.hh"
#include "info_impl_fw.hh"

namespace ost { namespace info {

class DLLEXPORT InfoGroup {
  friend class InfoHandle; // ctor access
  friend class InfoItem; // ctor access
public:

  //! retrieve parent group
  InfoGroup GetParent() const;

  //! set name
  void SetName(const String& name);
  //! retrieve name
  String GetName() const;

  //! retrieve path
  InfoPath GetPath() const;

  InfoGroup GetGroup(const InfoPath& path, bool use_defaults=true) const;
  
  InfoGroupList GetGroups(const InfoPath& path) const;
  
  InfoItemList GetItems(const InfoPath& path) const;
  
  InfoGroup CreateGroup(const String& name);
  bool HasGroup(const InfoPath& name, bool use_defaults=true) const;
  InfoGroup RetrieveGroup(const InfoPath& path, bool use_defaults=true);

  InfoItem GetItem(const InfoPath& path, bool use_defaults=true) const;
  InfoItem CreateItem(const String& name, const String& value);
  InfoItem CreateItem(const String& name, Real value);
  InfoItem CreateItem(const String& name, bool value);  
  InfoItem CreateItem(const String& name, int value);  
  InfoItem CreateItem(const String& name, const geom::Vec3& vector);
  bool HasItem(const InfoPath& path, bool use_defaults=true) const;
  InfoItem RetrieveItem(const InfoPath& path, bool use_defaults=true);

  void Remove(const InfoPath& path, bool use_defaults=false);
  void Remove(const InfoGroup& group);


  //! return attribute of given name
  String GetAttribute(const String& name) const;
  //! set attribute of given name, will be created if it does not yet exist
  void SetAttribute(const String& name, const String& value);
  //! returns true if attribute of this names exists
  bool HasAttribute(const String& name) const;
  //! remove attribute of given name
  void RemoveAttribute(const String& name);

  std::vector<String> GetAttributeList() const;

  String GetTextData() const;
  void SetTextData(const String& td);

  //! get comment just above item
  String GetComment() const;


  //! Apply visitor
  /*!
    descends through all items and subgroups
  */
  void Apply(InfoVisitor& v, bool visit_this=true);
  void Apply(InfoConstVisitor& v, bool visit_this=true) const;

  bool operator==(const InfoGroup& ref) const;
  bool operator!=(const InfoGroup& ref) const;
private:
  void Swap(InfoGroup& group);
  InfoGroup(const InfoHandle& root, const ElePtr& impl);

  bool do_group_lookup(std::vector<String>::const_iterator& pos,
                       const std::vector<String>::const_iterator& end,
                       InfoGroupList& subgroups) const;
  bool do_group_lookup(std::vector<String>::const_iterator& pos,
                       const std::vector<String>::const_iterator& end,
                       InfoGroup& subgroup) const;                       

  InfoGroup group_create(InfoGroup group,
                         std::vector<String>::const_iterator& pos,
                         const std::vector<String>::const_iterator& end);

  InfoHandle root_;
  ElePtr impl_;
};

}} // ns

#endif
