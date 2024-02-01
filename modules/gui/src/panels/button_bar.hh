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
#ifndef OST_GUI_BOTTOM_BAR_BUTTON_BAR
#define OST_GUI_BOTTOM_BAR_BUTTON_BAR


#include <ost/gui/module_config.hh>
#include <ost/gui/widget.hh>

#include "button_box.hh"
#include "event_button.hh"

#include <QWidget>
#include <QString>
#include <QHBoxLayout>
#include <QList>
namespace ost { namespace gui {

/// \brief button box
class DLLEXPORT_OST_GUI ButtonBar : public Widget {
  Q_OBJECT
public:
  ButtonBar(QWidget* parent=NULL);
  
  virtual bool Save(const QString& prefix);
  virtual bool Restore(const QString& prefix);

  virtual void AddButtonBox(ButtonBox* button_box_);

private slots:
  void ButtonDragged(EventButton* widget);
  void ButtonDropped(EventButton* widget);

private:
  QHBoxLayout* layout_;
  ButtonBox* drag_box_;
  EventButton* drag_button_;
};

}}

#endif
