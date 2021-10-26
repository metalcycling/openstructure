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
#ifndef OST_GUI_SCENE_WIN_CONTEXT_MENU_HH
#define OST_GUI_SCENE_WIN_CONTEXT_MENU_HH

/*
  Author: Stefan Scheuber
 */


#include <ost/gui/module_config.hh>
#include <ost/gui/scene_win/scene_win_model.hh>

#include <QFlags>
#include <QVector>
#include <QAbstractItemModel>
#include <QTreeView>
namespace ost { namespace gui {

enum ContextActionType
{
  GFX_NODE=0x1,
  GFX_OBJECT=0x2,
  ENTITY=0x4,
  ENTITY_VIEW=0x8,
  CUSTOM_VIEW=0x10,
  NOT_VISIBLE=0x20,
  NOT_HIDDEN=0x40,
  NOT_SCENE=0x80,
  SINGLE=0x100,
  MULTI=0x200,
  VIEWS_SAME_OBJECT=0x400,
  MAP=0x800,
  MAP_ORIGINAL=0x1000,
  MAP_DOWNSAMPLED=0x2000,
  MAP_DSAMPLED_AVAIL=0x4000
};

Q_DECLARE_FLAGS(ContextActionTypes, ContextActionType)


class DLLEXPORT_OST_GUI ContextMenu : public QObject {
  Q_OBJECT
public:
  ContextMenu(QTreeView* view, SceneWinModel* model);
  ~ContextMenu(){};

  void ShowMenu(const QPoint& pos);
  void Rename(QModelIndex index);

  void AddAction(QAction* action,ContextActionTypes type);

private slots:
  void AddView();
  void Rename();
  void DeleteView();

private:
  QVector<std::pair<QAction*, ContextActionTypes> > actions_;

  QTreeView* view_;
  SceneWinModel* model_;
};

}}

Q_DECLARE_OPERATORS_FOR_FLAGS(ost::gui::ContextActionTypes)

#endif
