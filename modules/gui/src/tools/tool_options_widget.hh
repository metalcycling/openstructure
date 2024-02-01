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
#ifndef OST_GUI_TOOL_OPTIONS_WIDGET_HH
#define OST_GUI_TOOL_OPTIONS_WIDGET_HH

/*
  Author: Marco Biasini
 */

#include <ost/gui/module_config.hh>

#include <ost/gui/tools/tool_options.hh>

#include <QWidget>

namespace ost { namespace gui {

class DLLEXPORT_OST_GUI EnumOptBinder : public QObject {
  Q_OBJECT
public:
  EnumOptBinder() {}
  EnumOptBinder(QObject* parent, ToolOptionEnum* enum_opt):
    QObject(parent), enum_opt_(enum_opt)
  { }
public slots:
  void IndexChanged(int index)
  {
    enum_opt_->SetIndex(index);
  }
private:
  ToolOptionEnum* enum_opt_;
};

class DLLEXPORT_OST_GUI IntOptBinder : public QObject {
  Q_OBJECT
public:
  IntOptBinder() {}  
  IntOptBinder(QObject* parent, ToolOptionInt* int_opt):
    QObject(parent), int_opt_(int_opt)
  { }    
public slots:
  void TextChanged(const QString& text)
  {
    int_opt_->SetValue(text.toInt());
  }
private:
  ToolOptionInt* int_opt_;
};

class DLLEXPORT_OST_GUI FloatOptBinder : public QObject {
  Q_OBJECT
public:
  FloatOptBinder() {}  
  FloatOptBinder(QObject* parent, ToolOptionFloat* float_opt):
    QObject(parent), float_opt_(float_opt)
  { }
public slots:
  void TextChanged(const QString& text)
  {
    float_opt_->SetValue(text.toFloat());
  }
private:
  ToolOptionFloat* float_opt_;
};
  
class DLLEXPORT_OST_GUI ToolOptionsWidget : public QWidget {
public:
  ToolOptionsWidget(ToolOptions* options, QWidget* parent=NULL);
private:
  void Populate();
  QWidget* MakeIntWidget(ToolOptionInt* opts);
  QWidget* MakeEnumWidget(ToolOptionEnum* opts);  
  QWidget* MakeFloatWidget(ToolOptionFloat* opts);  
  QWidget* MakeButtonWidget(ToolOptionButton* opts);
  ToolOptions* options_;
};

}}

#endif
