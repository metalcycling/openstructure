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
#ifndef OST_GUI_WIDGET_REGISTRY_HH
#define OST_GUI_WIDGET_REGISTRY_HH

#include <QWidget>
#include <QMap>

#include <ost/gui/module_config.hh>

/*
  Author: Marco Biasini
 */
namespace ost { namespace gui {

class Widget;

class DLLEXPORT_OST_GUI WidgetFactory : public QObject {
protected:
  WidgetFactory(const QString& id, const QString& name):
    id_(id), name_(name)
  {}
public:
  const QString& GetUniqueID() const { return id_; }
  const QString& GetFullName() const { return name_; }
  
  virtual Widget* Create(QWidget* parent)=0;
private:
  QString id_;
  QString name_;
};

#define OST_REGISTER_WIDGET_WITH_DEFAULT_FACTORY(ns, class_name, full_name)    \
  class class_name##Factory : public WidgetFactory {                           \
  public:                                                                      \
    class_name##Factory(): WidgetFactory(#ns"::"#class_name, full_name) {}     \
    Widget* Create(QWidget* parent)                                            \
    {                                                                          \
      return new class_name(parent);                                           \
    }                                                                          \
  };                                                                           \
  bool class_name##_registered=WidgetRegistry::Instance()->                    \
                               RegisterWidgetFactory(new class_name##Factory)
                               
#define OST_REGISTER_WIDGET(class_name, factory)                               \
  bool class_name##_registered=WidgetRegistry::Instance()->                    \
                               RegisterWidgetFactory(new factory)


/// \brief central registry for widgets
class DLLEXPORT_OST_GUI WidgetRegistry : public QObject {
public:
  static WidgetRegistry* Instance();
  /// \brief register widget factory
  /// 
  /// The widget registry takes ownership of the widget factory
  bool RegisterWidgetFactory(WidgetFactory* factory);

  /// \brief unregister widget factory
  ///
  /// The ownership of the widget factory is set to NULL
  bool UnregisterWidgetFactory(const QString& id);

  /// \brief create widget by id
  Widget* Create(const QString& id, QWidget* parent);
  QString GetFullName(const QString& id);
private:
  WidgetRegistry() {}
  QMap<QString, WidgetFactory*> id_to_factory_;
};
 

}}

#endif
