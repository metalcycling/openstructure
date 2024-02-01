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
#ifndef OST_SEQUENCE_VIEWER_SEQUENCE_TABLE_VIEW
#define OST_SEQUENCE_VIEWER_SEQUENCE_TABLE_VIEW

/*
  Author: Stefan Scheuber
 */


#include <ost/gui/module_config.hh>

#include "sequence_delegate.hh"

#include <QTableView>
#include <QWheelEvent>
namespace ost { namespace gui {

/// \brief QTableView with first column not moving
class DLLEXPORT_OST_GUI SequenceTableView : public QTableView {
  Q_OBJECT
public:
  SequenceTableView(QAbstractItemModel * model);
  QTableView* GetStaticRow();
  QTableView* GetStaticColumn();
  QTableView* GetStaticField();
  ~SequenceTableView();


  void SetSingleAlignment(bool flag) { single_alignment_=flag; }
  
  bool IsSingleAlignment() const { return single_alignment_; }
signals:
  void MouseWheelEvent(QWheelEvent* event);
  void CopyEvent(QKeyEvent* event);
  void AlignmentChanged();
public slots:
  void columnCountChanged(const QModelIndex& index, int old_count, int new_count);
  void rowCountChanged(const QModelIndex& index, int old_count, int new_count);
  void resizeColumnsToContents();
  void resizeRowsToContents();

protected:
  virtual void mouseDoubleClickEvent(QMouseEvent* event);
  virtual void mouseReleaseEvent(QMouseEvent* event);
  virtual void resizeEvent(QResizeEvent* event);
  virtual void wheelEvent(QWheelEvent* event);
  virtual void keyPressEvent(QKeyEvent* event);
  virtual QModelIndex moveCursor(CursorAction cursorAction, Qt::KeyboardModifiers modifiers);
  void scrollTo (const QModelIndex & index, ScrollHint hint = EnsureVisible);
  virtual int sizeHintForColumn(int column) const;
  virtual int sizeHintForRow(int row) const;

private slots:
  void InitStaticColumn();
  void InitStaticRow();
  void InitStaticField();
  void ResizeWidth(int index, int, int size);
  void ResizeHeight(int index, int, int size);

private:
  std::pair<int, int> GetSimpleSelection(int row);
  QTableView* static_column_;
  QTableView* static_row_;
  QTableView* static_field_;
  SequenceDelegate* delegate_;
  QModelIndex last_double_click_;
  bool        single_alignment_;
  void updateStaticColumn();
  void updateStaticRow();
  void updateStaticField();
};

}}

#endif
