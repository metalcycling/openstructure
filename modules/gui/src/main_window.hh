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
  Author: Ansgar Philippsen
*/

#ifndef OST_GUI_MAIN_WINDOW_HH
#define OST_GUI_MAIN_WINDOW_HH

#include <map>


#include <ost/base.hh>
#include <ost/gui/module_config.hh>


#include <QMainWindow>
#include <QMenu>

namespace ost { namespace gui {

class DLLEXPORT_OST_GUI MainWindow: public QMainWindow
{
  Q_OBJECT;

  typedef std::map<QWidget*,QDockWidget*> DockWidgetMap;
public:
  MainWindow(QWidget* parent=0, Qt::WindowFlags flags=Qt::WindowFlags());

  //! add a custom docking widget
  void AddDockWidget(QWidget* w, const QString& name, bool show=true, int area=2);
  //! remove a previously added custom docking widget
  void RemoveDockWidget(QWidget* w);

public slots:
  void OnDockMenu(QAction* a);
  void OnDockClose(QDockWidget* dock);

protected:
  QMenu* WindowMenu();
  const QMenu* WindowMenu() const;

private:
  QMenu* window_menu_;
  DockWidgetMap dock_widget_map_;
};

}} // ns

#endif
