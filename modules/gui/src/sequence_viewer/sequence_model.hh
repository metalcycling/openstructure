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
#ifndef OST_SEQUENCE_VIEWER_SEQUENCE_MODEL
#define OST_SEQUENCE_VIEWER_SEQUENCE_MODEL

/*
  Author: Stefan Scheuber
 */

#ifndef Q_MOC_RUN

#include <ost/mol/chain_view.hh>

#include <ost/seq/sequence_list.hh>
#include <ost/seq/alignment_handle.hh>

#include <ost/gfx/entity.hh>

#include "base_view_object.hh"
#include "alignment_view_object.hh"
#include "sequence_view_object.hh"

#include <QAbstractTableModel>
#include <QItemSelection>
#endif

namespace ost { namespace gui {

class SequenceModel : public QAbstractTableModel
{
  Q_OBJECT

public:
  SequenceModel(QObject *parent = 0);

  void InsertAlignment(const seq::AlignmentHandle& alignment);
  void InsertGfxEntity(const gfx::EntityP& entity);
  void InsertChain(QString& name, mol::ChainView& view);
  void InsertSequence(QString& name, seq::SequenceHandle& seq);
  void InsertSequences(const QList<QString>& names, seq::SequenceList& list);

  void RemoveAlignment(const seq::AlignmentHandle& alignment);
  void RemoveGfxEntity(const gfx::EntityP& entity);

  QModelIndexList GetModelIndexes(gfx::EntityP& entity, const mol::EntityView& view);
  QModelIndexList GetModelIndexes(const QString& subject, const QString& sequence_name=QString());

  int GetGlobalRow(BaseViewObject* obj, int row) const;

  const QStringList& GetDisplayModes();
  const QStringList& GetDisplayModes(const gfx::EntityP& entity);
  const QStringList& GetDisplayModes(const seq::AlignmentHandle& alignment);
  const QString& GetCurrentDisplayMode();
  const QString& GetCurrentDisplayMode(const gfx::EntityP& entity);
  const QString& GetCurrentDisplayMode(const seq::AlignmentHandle& alignment);
  void SetDisplayMode(const QString& mode);
  void SetDisplayMode(const gfx::EntityP& entity, const QString& mode);
  void SetDisplayMode(const seq::AlignmentHandle& alignment, const QString& mode);

  const PainterList& GetPainters(const QModelIndex& index) const;

  // abstract item model interface
  int rowCount(const QModelIndex& parent=QModelIndex()) const;

  int columnCount(const QModelIndex& parent=QModelIndex()) const;

  QVariant data(const QModelIndex& index, int role=Qt::DisplayRole) const;

  QVariant headerData(int section, Qt::Orientation orientation,
                      int role=Qt::DisplayRole) const;

  virtual Qt::ItemFlags flags(const QModelIndex& index=QModelIndex()) const;
  QPair<seq::AlignmentHandle, int> GetAlignmentForRow(int row);
  
  void EmitRowChanged(int row);
  
  void Clear();
public slots:
  void ZoomIn();
  void ZoomOut();
  void DoubleClicked(const QModelIndex& index);
  void SelectionChanged(const QItemSelection& sel, const QItemSelection& desel);
private:
  int GetColumnCount() const;
  SequenceViewObject* GetItem(const gfx::EntityP& entity);
  AlignmentViewObject* GetItem(const seq::AlignmentHandle& alignment);
  BaseViewObject* GetItem(const QModelIndex& index) const;
  QPair<int, BaseViewObject*> GetRowWithItem(int row) const;
  QPair<int, BaseViewObject*> GetRowWithItem(const QModelIndex& index) const;

  QList<BaseViewObject*> objects_;
  int max_columns;
  PainterList empty_painter_list_;
  QString empty_string_;
  QStringList empty_string_list_;
  QStringList display_modes_;
  QString current_display_mode_;
};


}}

#endif
