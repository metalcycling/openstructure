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
#ifndef OST_GUI_SCENE_SELECTION_HH
#define OST_GUI_SCENE_SELECTION_HH


#include <ost/mol/query_view_wrapper.hh>

#include <ost/gfx/entity_fw.hh>
#include <ost/gfx/gfx_node_fw.hh>
#include <ost/gui/module_config.hh>

#include <QObject>
#include <QMenu>
#include <QPoint>
namespace ost { namespace gui {

class DLLEXPORT_OST_GUI SceneSelection : public QObject {
  Q_OBJECT
public:
  static SceneSelection* Instance();
  int GetActiveNodeCount() const;
  gfx::GfxNodeP GetActiveNode(unsigned int pos) const;
  int GetActiveViewCount() const;
  mol::EntityView GetActiveView(unsigned int pos) const;
  gfx::EntityP GetViewEntity() const;

public slots:
  void CenterOnObjects();
  void Delete();
  void ViewDensitySlices();
  void ShowDownsampledMap();
  void ShowOriginalMap();  
  void Deselect();
  void Show();
  void Hide();
  void SelectAllViews();
  void DeselectAllViews();
  void MakeVisible();
  void MakeHidden();
  void ShowExclusive();
  void HideExclusive();
  mol::EntityView GetViewUnion();


private slots:
  void SetActiveNodes(gfx::NodePtrList nodes, gfx::EntityP entity, mol::QueryViewWrapperList views);

private:
  SceneSelection();
  gfx::NodePtrList nodes_;
  gfx::EntityP view_entity_;
  mol::QueryViewWrapperList views_;
  static SceneSelection* scene_selection_;
  QMenu* menu_;
};

}}

#endif
