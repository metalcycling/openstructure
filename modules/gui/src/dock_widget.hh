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

#ifndef OST_GUI_DOCK_WIDGET_HH
#define OST_GUI_DOCK_WIDGET_HH


#include <ost/gui/module_config.hh>

#include <QDockWidget>
namespace ost { namespace gui {

class DLLEXPORT_OST_GUI DockWidget: public QDockWidget
{
  Q_OBJECT;
public:
  DockWidget(const QString& title, QWidget* parent=0, Qt::WindowFlags flags=Qt::WindowFlags());
  virtual void closeEvent(QCloseEvent* e);
signals:
  void OnClose(QDockWidget* dock);
};

}} // ns

#endif
