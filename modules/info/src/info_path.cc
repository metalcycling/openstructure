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
  info path

  Author: Ansgar Philippsen
*/

#include <iostream>
#include <boost/tokenizer.hpp>

#include "info_path.hh"

namespace ost { namespace info {

InfoPath::InfoPath(const String& path)
{
  init(path);
}

InfoPath::InfoPath(const char* cp)
{
  init(cp);
}

InfoPath::InfoPath(const std::vector<String> lst,bool rf):
  list_(lst),
  relative_flag_(rf)
{}

void InfoPath::init(const String& path)
{
  list_.clear();
  if(path.empty()) {
    relative_flag_=false;
    return;
  }
  relative_flag_=true;
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  static boost::char_separator<char> sep("/","",boost::keep_empty_tokens);
  tokenizer tokens(path,sep);
  tokenizer::iterator it = tokens.begin();
  if(it->empty()) {
    // root
    relative_flag_=false;
    ++it;
  }
  for (;it!=tokens.end(); ++it) {
    if(!it->empty()) {
      list_.push_back(*it);
    }
  }
}

std::ostream& operator<<(std::ostream& os, const InfoPath& p)
{
  std::vector<String> plist=p.GetList();
  if(!p.IsRelative()) {
    os << "/";
  }
  for(std::vector<String>::const_iterator it=plist.begin();it!=plist.end();++it) {
    os << *it << "/";
  }
  
  return os;
}

InfoPath InfoPath::Strip() const
{
  if(list_.size()>1) {
    std::vector<String> plist(list_.begin(),list_.end()-1);
    return InfoPath(plist,relative_flag_);
  } else {
    return InfoPath("");
  }
}

InfoPath operator+(const InfoPath& p1, const InfoPath& p2)
{
  std::vector<String> pathlist(p1.list_.begin(),p1.list_.end());
  pathlist.insert(pathlist.end(),p2.list_.begin(),p2.list_.end());
  return InfoPath(pathlist,p1.IsRelative());
}

}} // ns
