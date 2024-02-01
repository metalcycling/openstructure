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
#ifndef OST_SEQUENCE_VIEWER_BASE_VIEW_OBJECT
#define OST_SEQUENCE_VIEWER_BASE_VIEW_OBJECT

/*
  Author: Stefan Scheuber
 */


#include "base_row.hh"

#include <QObject>
#include <QList>

namespace ost { namespace gui {

class BaseViewObject : public QObject
{
  Q_OBJECT


public:
  BaseViewObject(QObject* parent = 0);

  void InsertRow(int pos, BaseRow* row);
  void RemoveRow(BaseRow* row);
  BaseRow* GetRow(int pos);
  int GetRowCount();
  int GetMaxColumnCount() const;

  virtual void SetSelection(int row, const QSet<int>& added, const QSet<int>& removed);

  virtual QVariant GetData(int row, int column, int role);
  virtual bool SetData(int row, int column, const QVariant& value, int role);
  virtual Qt::ItemFlags Flags(int row, int column) const;

  virtual const QStringList& GetDisplayModes();
  virtual const QString& GetCurrentDisplayMode();
  virtual void SetDisplayMode(const QString& mode);

  void DoubleClicked(int row, int column);
  void ZoomIn();
  void ZoomOut();

protected:
  virtual void AddDisplayMode(const QString& mode);

  QList<BaseRow*> rows_;
  QString current_display_mode_;
  QStringList display_modes_;
};


}}

#endif
