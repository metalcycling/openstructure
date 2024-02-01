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
#ifndef OST_SEQUENCE_VIEWER_ALIGNMENT_VIEW_OBJECT
#define OST_SEQUENCE_VIEWER_ALIGNMENT_VIEW_OBJECT

/*
  Author: Stefan Scheuber
 */

#ifndef Q_MOC_RUN

#include <ost/seq/alignment_handle.hh>

#include <ost/gfx/gradient.hh>

#include "sequence_view_object.hh"
#endif

namespace ost { namespace gui {

class AlignmentViewObject : public SequenceViewObject
{
  Q_OBJECT

public:
  AlignmentViewObject(const seq::AlignmentHandle& alignment, QObject* parent = 0);

  QVariant GetData(int row, int column, int role);
  const seq::AlignmentHandle& GetAlignment();

  void SetDisplayMode(const QString& mode);

private:
  gfx::Gradient gradient_;
  seq::AlignmentHandle alignment_;
  QMap<int, QColor> conservation_1_;
  QMap<int, QColor> conservation_2_;

  static QMap<QString,int> group_map_;

  static const QString conservation_mode_1;
  static const QString conservation_mode_2;
};


}}

#endif
