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
#ifndef OST_SEQUENCE_VIEWER_TITLE_ROW
#define OST_SEQUENCE_VIEWER_TITLE_ROW

/*
  Author: Stefan Scheuber
 */


#include "base_row.hh"
#include <QObject>

namespace ost { namespace gui {

class TitleRow : public BaseRow
{
  Q_OBJECT

public:
  TitleRow(QObject* parent);

  virtual QVariant GetData(int column, int role) const;
  virtual Qt::ItemFlags Flags(int column) const;
  virtual void DoubleClicked(int column);
private:
  QFont font_;
  QFont small_font_;
};
}}

#endif
