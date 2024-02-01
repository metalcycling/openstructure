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
#ifndef OST_GUI_WIDGET_POOL_HH
#define OST_GUI_WIDGET_POOL_HH

#include <ost/gui/widget.hh>
#include <QSettings>
/*
  Author: Marco Biasini, Stefan Scheuber
 */
namespace ost { namespace gui {

struct PooledWidget
{
  PooledWidget():
    name(""), available(-1)
  { }

  PooledWidget(const QString& n, int lim=-1):
    name(n), available(lim), limit(lim)
  { }
  QString name;
  int     available;
  int     limit;
};

/// \brief This Widget pool manages instances of widgets
/// Widgets can be added to this pool by simply passing the classname and the amount of available instances. This class is coupled to the WidgetRegistry .
class DLLEXPORT_OST_GUI WidgetPool : public QObject {
  Q_OBJECT
public:
  WidgetPool(QObject* parent=NULL);

  /// \brief add a Widget to the WidgetPool
  /// This method adds a class to this widget pool. It does NOT add the given class to the WidgetRegistry.
  /// You should first add the class to the WidgetRegistry use the \ref WidgetRegistry::RegisterWidgetFactory(WidgetFactory* factory) for this.
  /// \param name full qualified class name. Is used to identify the class.
  /// \param limit amount of available instances of the class in the WidgetRegistry. If unlimited then this can be set to -1.
  void Add(const QString& name, int limit=-1);

  /// \brief returns whether this class has available instances or not
  /// \param name Full qualified class of the class you would like to check
  /// \return returns true if there are available instances. False if all instances are used or class is not in widget registry.
  bool IsAvailable(const QString& name);

  /// \brief returns a instance of the given class name
  /// It does not \ref WidgetPool::Take(const QString& name) a instance out of the pool. The user of this class is responsible to call the take method before calling this method.
  /// \param name Full qualified class of the class of which you would like to have an instance.
  /// \return Returns a instance of the desired class. Returns NULL if all instances are used or class is not in widget registry.
  Widget* Create(const QString& name);

  /// \brief Takes a instance of a class out of the pool
  /// This method takes one available instances out of the pool. If there are no available instances or unlimited amount of instances, nothing will be changed. Call this method always before creating an instance. If false is returned do not call \ref WidgetPool::Create(const QString& name) !
  /// \param name full qualified class of the class of which you would like to take an instance.
  /// \return returns false if no instance is available or class is not known. Otherwise true will be returned.
  bool Take(const QString& name);

  /// \brief Returns a instance of a class into the pool
  /// This method returns one instance back into the pool. If there is an unlimited amount of instances of this class, nothing will be changed inside the pool.
  /// \param name full qualified class of the class of which you give back an instance.
  /// \return returns false if class is not known. Otherwise true will be returned.
  bool Give(const QString& name);

  /// \brief Get a List with all available Widgets
  /// \return The returned list contains all widgets which have available instances in the pool. Instances which have no available instances won't be in this list.
  QList<QString> GetAvailableWidgets();

private:
  QMap<QString, PooledWidget> pool_;
};

}}
#endif /* OST_GUI_WIDGET_POOL_HH */
