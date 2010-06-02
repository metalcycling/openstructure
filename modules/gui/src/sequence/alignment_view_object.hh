//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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
#ifndef OST_SEQUENCE_VIEWER_ALIGNMENT_VIEW_OBJECT
#define OST_SEQUENCE_VIEWER_ALIGNMENT_VIEW_OBJECT

/*
  Author: Stefan Scheuber
 */

#include <ost/seq/alignment_handle.hh>

#include "sequence_view_object.hh"

namespace ost { namespace gui {

class AlignmentViewObject : public SequenceViewObject
{
  Q_OBJECT

public:
  AlignmentViewObject(const seq::AlignmentHandle& alignment, QObject* parent = 0);

  QVariant GetData(int row, int column, int role);

private:
  seq::AlignmentHandle alignment_;
  QMap<int, QColor> conservation_;

  static QMap<QString,int> group_map_;
};


}}

#endif
