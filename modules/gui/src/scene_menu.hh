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
#ifndef OST_GUI_SCENE_MENU_HH
#define OST_GUI_SCENE_MENU_HH


#include <ost/mol/atom_handle.hh>
#include <ost/gfx/gfx_object_fw.hh>

#include <QObject>
#include <QPoint>
#include <QMenu>
#include <QString>
namespace ost { namespace gui {

class SceneMenu : public QObject{
  Q_OBJECT
public:
  SceneMenu();
  virtual ~SceneMenu();
  void SetPickPoint(const QPoint& pos);

public slots:
  void ShowMenu(const QPoint& pos);
  void SelectBackbone();
  void SelectSidechains();
  void SelectLigands();
  void SelectAroundAtom();
  void SelectAroundResidue();
  void SelectAtom();
  void SelectResidue();
  void Select(const String& str);
  void UnSelect();
  void CenterOnAtom();
  void CenterOnResidue();
  void CenterOnChain();
  void CenterOnView();
  void CenterOnSelection();

private:
  void PrepareContextMenu(const QPoint& pos);
  std::pair<gfx::GfxObjP, mol::AtomHandle>   picked_;
  QMenu* context_menu_;
  QMenu* select_menu_;
  QMenu* center_menu_;
};

}} // ns
#endif //OST_GUI_SCENE_MENU_HH
