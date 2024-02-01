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
#ifndef OST_SEQUENCE_VIEWER_BASE_ROW
#define OST_SEQUENCE_VIEWER_BASE_ROW

/*
  Author: Stefan Scheuber
 */


#include "painter.hh"

#include <QObject>
#include <QModelIndex>
#include <QList>
namespace ost { namespace gui {

class BaseRow : public QObject
{
  Q_OBJECT

public:
  BaseRow(QObject *parent = 0);
  BaseRow(QFont font, QObject *parent = 0);

  void Init();

  virtual int GetColumnCount() const;

  void InsertPainter(Painter* painter, int pos = -1);
  void RemovePainter(Painter* painter);
  Painter* GetPainter(int pos);
  int GetPainterCount();
  const PainterList& GetPainters() const;

  const QFont& GetFont() const;
  void SetFont(const QFont& font);
  const QSize& GetFontSize() const;

  virtual const QSize& GetCellSize() const;

  virtual QVariant GetData(int column, int role) const;
  virtual bool SetData(int column, const QVariant& value, int role);
  virtual Qt::ItemFlags Flags(int column) const;
  virtual void DoubleClicked(int column);
  virtual void SetSelection(const QSet<int>& added, const QSet<int>& removed);

  virtual void ZoomIn();
  virtual void ZoomOut();

private:
  bool IsPainterPosValid(int pos);
  PainterList painter_;
  QFont font_;
  QSize default_font_size_;
  QSize default_cell_size_;
};

typedef QList<BaseRow*> BaseRowList;

}}

#endif
