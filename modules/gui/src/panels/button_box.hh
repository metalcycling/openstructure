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
#ifndef OST_GUI_BOTTOM_BAR_BUTTON_BOX
#define OST_GUI_BOTTOM_BAR_BUTTON_BOX


#include <ost/gui/module_config.hh>
#include <ost/gui/widget.hh>

#include "event_button.hh"

#include <QWidget>
#include <QString>
#include <QGridLayout>
#include <QList>
#include <QLabel>
namespace ost { namespace gui {

/// \brief button box
class DLLEXPORT_OST_GUI ButtonBox : public Widget {
  Q_OBJECT
public:
  ButtonBox(const QString& label, QWidget* parent=NULL);
  
  virtual bool Save(const QString& prefix);
  virtual bool Restore(const QString& prefix);
  virtual void AddButton(const QString& name, Widget* widget, bool pressed);
  virtual void InsertButton(EventButton* position, EventButton* new_button);
  virtual void RemoveButton(EventButton* button);
  virtual void MoveButton(EventButton* position, EventButton* mov_button);

  virtual void dropEvent(QDropEvent* event);
  virtual void dragEnterEvent (QDragEnterEvent* event );

signals:
  void ToggleWidget(Widget* widget);
  void Removed(Widget* widget);
  void Inserted(Widget* widget, int pos, bool isHidden);
  void ButtonDragged(EventButton* button);
  void ButtonDropped(EventButton* button);

private slots:
  void ButtonPressed(Widget* widget);

private:
  QHBoxLayout* layout_;
  QHBoxLayout* buttons_layout_;
  QLabel* label_;
  QList<EventButton *>* buttons_;
};

}}

#endif
