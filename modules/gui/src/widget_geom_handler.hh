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
#ifndef OST_GUI_WIDGET_GEOM_HANDLER_HH
#define OST_GUI_WIDGET_GEOM_HANDLER_HH


#include <ost/gui/module_config.hh>

#include <QWidget>
namespace ost { namespace gui {

class DLLEXPORT_OST_GUI WidgetGeomHandler: public QObject{
  Q_OBJECT;
  public:
    WidgetGeomHandler(const QString& name, QWidget* widget, QWidget* parent=NULL);

    void init(const QString& name, QWidget* widget);

    void LoadGeom(const QString& prefix);
    void SaveGeom(const QString& prefix);
    QWidget* GetWidget() { return widget_; }
  private:
    QWidget* widget_;
};
}}

#endif
