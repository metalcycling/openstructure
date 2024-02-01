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

#include "thin_splitter.hh"

#include <QRegion>
#include <QPainter>
#include <QResizeEvent>
#include <QPaintEvent>
#include <QSplitterHandle>
#include <QApplication>
#include <QPalette>
/*
   Author: Marco Biasini
 */


namespace ost { namespace gui {

namespace {

class ThinSplitterHandle : public QSplitterHandle {
public:
  ThinSplitterHandle(Qt::Orientation orientation, QSplitter* parent):
    QSplitterHandle(orientation, parent)
  {
    this->setMask(QRegion(contentsRect()));
    this->setAttribute(Qt::WA_MouseNoMask, true);
  }
  
  virtual void resizeEvent(QResizeEvent *event)
  {
    if (this->orientation()==Qt::Horizontal) {
#if defined(__APPLE__)
      this->setContentsMargins(1, 0, 0, 0);
#else
      this->setContentsMargins(1, 0, 1, 0);
#endif
    } else {
#if defined(__APPLE__)
      this->setContentsMargins(0, 1, 0, 0);
#else                          
      this->setContentsMargins(0, 1, 0, 1);
#endif
    }
    this->setMask(QRegion());
  }
  
  virtual void paintEvent(QPaintEvent* event)
  {
    QPainter p(this);
    QPalette palette=qApp->palette();
    p.setBrush(QBrush(QColor(120,120,120)));
    p.setPen(QColor(120,120,120));
    p.drawRect(this->contentsRect());
  }
};

}

ThinSplitter::ThinSplitter(Qt::Orientation orientation, QWidget* parent):
  QSplitter(orientation, parent)
{
  this->setChildrenCollapsible(false);    
  this->setHandleWidth(1);  
}

QSplitterHandle* ThinSplitter::createHandle()
{
  return new ThinSplitterHandle(this->orientation(), this);
}

}}
