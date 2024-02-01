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
#ifndef OST_GUI_PANEL_BAR_WIDGET_HOLDER_HH
#define OST_GUI_PANEL_BAR_WIDGET_HOLDER_HH


#include <ost/gui/module_config.hh>

#include "panel_bar.hh"

#include <QWidget>
#include <QMouseEvent>
class QToolBar;
class QAction;

/*
  Author: Marco Biasini
 */

namespace ost { namespace gui {

class SideBar;
class Widget;

class DLLEXPORT_OST_GUI PanelBarWidgetHolder : public QWidget {
  Q_OBJECT
public:
  PanelBarWidgetHolder(Widget* widget, Qt::Orientation orientation, QWidget* parent=NULL);
  ~PanelBarWidgetHolder();

  void mousePressEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);

public slots:
  void Remove();
public:
  Widget* GetWidget() { return widget_; }
  void SetWidget(Widget* widget);
private slots:  
  void LabelClick(bool checked);  
private:
  void SetupToolBar();
  QToolBar*      toolbar_;  
  Widget*        widget_;
  QAction*       label_;
  QPoint         drag_start_pos_;
};

}}


#endif
