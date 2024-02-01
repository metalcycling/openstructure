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
  info visitor interface
  
  Author: Ansgar Philippsen
*/

#ifndef OST_INFO_VISITOR_HH
#define OST_INFO_VISITOR_HH

#include <ost/base.hh>
#include "module_config.hh"
#include "info_fw.hh"

namespace ost { namespace info {

class DLLEXPORT_OST_INFO InfoVisitor {
public:
  virtual ~InfoVisitor();

  //! visit routine called for each item
  /*!
    empty per default
  */
  virtual void VisitItem(InfoItem& item);

  //! visit routine called for each subgroup 
  /*!
    return value indicates wether subgroups of this subgroup
    shall be visited - upon true (the default), the descent will
    continue recursively
  */
  virtual bool VisitGroup(InfoGroup& group);

  //! called after all subgroups have been visited
  virtual void VisitGroupFinish(InfoGroup& group);
};

class DLLEXPORT_OST_INFO InfoConstVisitor {
public:
  virtual ~InfoConstVisitor();

  //! visit routine called for each item
  /*!
    empty per default
  */
  virtual void VisitItem(const InfoItem& item);

  //! visit routine called for each subgroup 
  /*!
    return value indicates wether subgroups of this subgroup
    shall be visited - upon true (the default), the descent will
    continue recursively
  */
  virtual bool VisitGroup(const InfoGroup& group);

  //! called after all subgroups have been visited
  virtual void VisitGroupFinish(const InfoGroup& group);
};

}} // ns


#endif
