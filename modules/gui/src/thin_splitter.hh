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
#ifndef OST_GUI_THIN_SPLITTER_HH
#define OST_GUI_THIN_SPLITTER_HH

#include <QSplitter>

/*
  Author: Marco Biasini
 */

class QSplitterHandle;

namespace ost { namespace gui {

/// \brief splitter with 1 pixel wide handles.
class ThinSplitter : public QSplitter {
public:
  ThinSplitter(Qt::Orientation orientation, QWidget* parent=0);
protected:
  QSplitterHandle* createHandle();
};
 
}}


#endif
