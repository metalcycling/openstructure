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
  verbose info visitor, a test implementation of the info visitor
  
  Author: Ansgar Philippsen
*/

#ifndef OST_INFO_VERBOSE_VISITOR_HH
#define OST_INFO_VERBOSE_VISITOR_HH

#include "info_visitor.hh"

namespace ost { namespace info {

class DLLEXPORT_OST_INFO VerboseInfoVisitor: public InfoVisitor {
public:
  virtual void VisitItem(InfoItem& item);
  virtual bool VisitGroup(InfoGroup& group);
};

}} // ns


#endif
