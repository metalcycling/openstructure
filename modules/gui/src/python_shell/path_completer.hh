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
#ifndef PATH_COMPLETER_HH
#define PATH_COMPLETER_HH

/*
  Authors: Marco Biasini, Andreas Schenk
 */

#include "completer_base.hh"


namespace ost { namespace gui {

//fw decl
class OstDirModel;

class PathCompleter : public OstCompleterBase
{
Q_OBJECT
public:
  PathCompleter(QObject * parent = 0);
  ~PathCompleter();
protected:
  OstDirModel* model_;
};

}}//ns

#endif
