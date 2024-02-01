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
#include "diag.hh"

namespace ost { namespace conop {

String Diag::Format(bool colored) const
{
  std::stringstream ss;
  std::vector<String> strings;
  for (std::vector<ArgDesc>::const_iterator
       i=args_.begin(), e=args_.end(); i!=e; ++i) {
   switch (i->type) {
     case DIAG_ARG_TYPE_ATOM:
       if (colored) {
         strings.push_back("\033[0;32m"+
                           atoms_[i->index].GetQualifiedName()+"\033[0m");
       } else {
         strings.push_back(atoms_[i->index].GetQualifiedName());
       }
       break;
     case DIAG_ARG_TYPE_RESIDUE:
        if (colored) {  
          strings.push_back("\033[0;32m"+
                            residues_[i->index].GetQualifiedName()+"\033[0m");
       } else {
          strings.push_back(residues_[i->index].GetQualifiedName());
       }
       break;
     case DIAG_ARG_TYPE_CHAIN:
        if (colored) {    
          strings.push_back("\033[0;32m"+
                            chains_[i->index].GetName()+"\033[0m");
       } else {
         strings.push_back(chains_[i->index].GetName());
       }
       break;
     case DIAG_ARG_TYPE_STRING:
       strings.push_back(strings_[i->index]);
       break;
     case DIAG_ARG_TYPE_INT:
       ss << ints_[i->index];
       strings.push_back(ss.str());
       ss.str("");
       break;
   }
  }
  for (size_t i=0; i<format_.size(); ++i) {
    if (format_[i]=='%') {
      // scan for number
      ++i;
      assert(i<format_.size());
      if (format_[i]=='%' ) {
        ss << '%';
        continue;
      }
      bool plural_s=false;
      if (format_[i]=='s') {
        plural_s=true;
        ++i;
      }
      const char *start=&format_[i];
      char *end=NULL;
      long int id=strtol(start, &end, 10);
      assert(start!=end);
      assert(id>=0 && static_cast<unsigned long>(id)<strings.size());
      if (plural_s) {
        if (ints_[args_[id].index]!=1) {
          ss << "s";
        }
      } else {
        ss << strings[id];
      }
      i+=end-start-1;
      continue;
    } 
    ss << format_[i];
  }
  return ss.str();
}

}} /* ost::conop */
