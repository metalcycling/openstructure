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
#ifndef OST_GUI_SCENE_WIN_SCENE_WIN_HH
#define OST_GUI_SCENE_WIN_SCENE_WIN_HH


#include <ost/mol/query_view_wrapper.hh>

#include <ost/gfx/gfx_node_fw.hh>
#include <ost/gfx/entity_fw.hh>

#include <ost/gui/module_config.hh>
#include <ost/gui/widget.hh>

#include <QTreeView>
#include <QItemSelection>


/*
  Authors: Marco Biasini, Ansgar Philippsen, Stefan Scheuber
*/

namespace ost { namespace gui {

class SceneWinModel;
class ContextMenu;

// the display window for all graphical objects
class DLLEXPORT_OST_GUI SceneWin: public Widget
{
  Q_OBJECT;
public:
  SceneWin(QWidget* parent=NULL);
  ~SceneWin();

signals:
  void ActiveNodesChanged(gfx::NodePtrList nodes, gfx::EntityP entity, mol::QueryViewWrapperList views);
  void ActiveNodesChanged();

public:
  virtual bool Save(const QString& prefix) { return true; }
  virtual bool Restore(const QString& prefix) { return true; }

  void AddView(gfx::EntityP entity, mol::EntityView view);

  ContextMenu* GetContextMenu();

public slots:
  void OnSelectionChange(const QItemSelection& sel, const QItemSelection& desel);

  void ContextMenuRequested(const QPoint& pos);

  SceneWinModel* GetModel() const;

  void Update();

private slots:
  void RowsInserted(const QModelIndex & parent, int start, int end);
  void DoubleClicked(const QModelIndex & index);

private:
  SceneWinModel* model_;
  QTreeView* view_;
  ContextMenu* context_menu_;
};

}} // ns

#endif
