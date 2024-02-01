//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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

/*
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMG_GUI_DATAVIEWER_INFO_PANEL_H_
#define IMG_GUI_DATAVIEWER_INFO_PANEL_H_


#include <ost/img/data_types.hh>
#include <ost/img/extent.hh>
#include <ost/img/image_fw.hh>

#include <ost/gui/module_config.hh>

#include <QWidget>
#include <QLabel>
#include <QMenu>
namespace ost { namespace img { namespace gui {

class DLLEXPORT_OST_GUI InfoPanelLabel: public QLabel
{
  Q_OBJECT;
public:
  InfoPanelLabel(const QString& t="", QWidget* parent=0);
};

class InfoPanel: public QWidget
{
  Q_OBJECT;
public:
  InfoPanel(QWidget* parent);
  void SetImageInfo(const ConstData& data);
  void SetMousePoint(const geom::Vec2& point, Real val);
  void SetClickPoint(const geom::Vec2& point, Real val);
  void SetMousePoint(const geom::Vec2& point, Complex val);
  void SetClickPoint(const geom::Vec2& point, Complex val);
  void ClearClickPosition();

public slots:
  void SetSelection(const Extent& e);
  void ClearSelection(); 
private:
  void UpdateMouseAndDistances(const geom::Vec2& point);
  void UpdateClick(const geom::Vec2& point);
  InfoPanelLabel* image_extent_;
  InfoPanelLabel* image_size_;
  InfoPanelLabel* image_sampling_;
  InfoPanelLabel* image_type_;
  InfoPanelLabel* mouse_xy_;
  InfoPanelLabel* mouse_val_;
  InfoPanelLabel* click_xy_;
  InfoPanelLabel* click_val_;
  InfoPanelLabel* selection_xy_;
  InfoPanelLabel* size_xy_;
  InfoPanelLabel* distance_;
  QMenu* popupmenu_;
  bool amp_pha_;
  geom::Vec2 click_;
};

}}}  //ns

#endif
