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
#ifndef OST_GUI_WIDGET_HH
#define OST_GUI_WIDGET_HH


#include <ost/gui/module_config.hh>

#include <QWidget>
#include <QMetaType>
/*
  Author: Marco Biasini
 */

namespace ost { namespace gui {

typedef QList<QAction*> ActionList;

/// \brief Base class for widgets
/// 
/// The base class defines an interface for restoring and saving the state
/// of the internal widget. 
/// 
class DLLEXPORT_OST_GUI Widget : public QWidget {
  Q_OBJECT

public:
  Widget(QWidget* widget, QWidget* parent=NULL);
  
  /// \brief set whether the widget should be destroyed when closed
  void SetDestroyOnClose(bool flag);
  
  /// \brief whether the widget should be destroyed upon receiving a 
  ///     close event
  bool DestroyOnClose() const;
  
  /// \brief get Widget specific actions
  virtual ActionList GetActions();

  QWidget* GetInternalWidget();
  void SetInternalWidget(QWidget* widget);
  
  
  const QString& GetUniqueID() const
  {
    if (unique_id_.size()==0) {
      unique_id_=this->metaObject()->className();
    }
    return unique_id_;
  }
  
  void SetUniqueID(const QString& id)
  {
    unique_id_=id;
  }
  
  virtual bool Restore(const QString& prefix)=0;
  
  virtual bool Save(const QString& prefix)=0;
  
private:
  mutable QString  unique_id_;
  QWidget* internal_;
  bool     destroy_on_close_;
};

}}

Q_DECLARE_METATYPE(ost::gui::Widget*)

#endif
