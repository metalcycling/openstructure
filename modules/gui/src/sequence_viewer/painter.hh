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
#ifndef OST_SEQUENCE_VIEWER_PAINTER
#define OST_SEQUENCE_VIEWER_PAINTER

/*
  Author: Stefan Scheuber
 */


#include <ost/mol/alg/sec_structure_segments.hh>

#include <QObject>
#include <QPainter>
#include <QStyleOptionViewItem>
#include <QModelIndex>
#include <QVarLengthArray>
#include <QList>
namespace ost { namespace gui {

class Painter : public QObject
{
  Q_OBJECT

public:
  Painter(QObject* parent = 0):QObject(parent){};
  virtual void Paint(QPainter* painter, const QStyleOptionViewItem& option,
                           const QModelIndex& index) = 0;

};

typedef QList<Painter*> PainterList;

}}

Q_DECLARE_METATYPE(QVarLengthArray<ost::mol::SecStructure>)

#endif
