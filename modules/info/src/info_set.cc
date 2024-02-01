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
  Author: Ansgar Philippsen
*/

#include <iostream>

#include <boost/regex.hpp>

#include "info_set.hh"
#include "info_visitor.hh"
#include "info_group.hh"
#include "info_item.hh"

namespace ost { namespace info {

namespace {

// helper class, actual work horse for set extraction

class SetExtractor : public InfoConstVisitor
{
public:
  SetExtractor(): 
    init_(false),
    current_set_("")
  {
    //std::cerr << "entering set extractor" << std::endl;
  }

  virtual void VisitItem(const InfoItem& item) {
    if(init_ && item.GetName()=="entry") {
      //std::cerr << " " << item.GetValue() << std::endl;
      current_set_.Add(item.GetValue());
    }
  }

  virtual bool VisitGroup(const InfoGroup& group) {
    if(group.GetName()=="set") {
      //std::cerr << "set " << group.GetAttribute("name") << std::endl;
      init_=true;
      current_set_=InfoSet(group.GetAttribute("name"));
      return true; // descend into subgroups
    }
    return false; // dont descend into subgroups
  }

  virtual void VisitGroupFinish(const InfoGroup& group) {
    if (init_) {
      //std::cerr << "finished with set " << current_set_.GetName() << std::endl;
      set_map_[current_set_.GetName()]=current_set_;
      init_=false;
    }
  }
    
  InfoSetMap GetSets() {
    return set_map_;
  }

  bool init_;
  InfoSet current_set_;
  InfoSetMap set_map_;
};

}  // anon ns


// main extraction routine

InfoSetMap ExtractSets(const InfoGroup& ig)
{
  SetExtractor se;
  ig.Apply(se,false);
  return se.GetSets();
}


// info set class

InfoSet::InfoSet(const String& name):
  name_(name)
{}

String InfoSet::GetName() const {return name_;}

bool InfoSet::Match(const String& s) const
{
  for(std::vector<String>::const_iterator it=entry_list_.begin();
      it!=entry_list_.end();
      ++it) {
    boost::regex expr(*it);
    if(boost::regex_match(s,expr)) {
      return true;
    }
  }
  return false;
}

void InfoSet::Add(const String& entry) {entry_list_.push_back(entry);}

void InfoSet::Clear() {entry_list_.clear();}


}} // ns

