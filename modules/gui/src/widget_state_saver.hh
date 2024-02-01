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
#ifndef WIDGET_STATE_SAVER_HH_
#define WIDGET_STATE_SAVER_HH_

/*
   Author: Andreas Schenk
 */

#include <ost/gui/module_config.hh>

#include <QSettings>
#include <QPoint>
#include <QSize>
#include <QMainWindow>
#include <QDebug>

namespace ost{namespace gui{

class DLLEXPORT_OST_GUI WidgetStateSaverConfig
{
public:
  static void IgnoreSettings(bool flag=true);
  static bool SettingsIgnored();
protected:  
  static bool flag_; 
};

struct unused;
template<class BASE, typename T1=unused, typename T2=unused, typename T3=unused>
class WidgetStateSaver: public BASE
{
public:
  WidgetStateSaver(const QString& name, QWidget * parent ,T1 t1, T2 t2,T3 t3):
  BASE(parent,t1,t2,t3),
  name_(name),
  pos_(),
  size_(),
  visible_(true),
  default_(false)
  {
    set_state_();
  }
  
  WidgetStateSaver(const QString& name,QPoint pos, QSize size, QWidget * parent ,T1 t1, T2 t2,T3 t3):
  BASE(parent,t1,t2,t3),
  name_(name),
  pos_(pos),
  size_(size),
  visible_(true),
  default_(true)
  {
    set_state_();
  }

  void ResetStateSaver(const QString& name)
  {
    name_=name;
    set_state_();
  }

protected:
  
  void set_state_()
  {
    QSettings settings;
    settings.beginGroup(name_);
    if(settings.contains("pos") && ! WidgetStateSaverConfig::SettingsIgnored()){
      BASE::move(settings.value("pos").toPoint());
      BASE::resize(settings.value("size").toSize());
      //BASE::setVisible(settings.value("visible").toBool());
      if(dynamic_cast<QMainWindow* >(this)){
        dynamic_cast<QMainWindow* >(this)->restoreState(settings.value("state").toByteArray());
      }
    }else if(default_){
      BASE::move(pos_);
      BASE::resize(size_);
      //BASE::setVisible(visible_);
    }
    settings.endGroup();
  }
  
  virtual ~WidgetStateSaver()
  {
    QSettings settings;
    settings.beginGroup(name_);
    settings.setValue("pos", BASE::pos());
    settings.setValue("size", BASE::size());
    //settings.setValue("visible", /*BASE::isVisible()*/visible_);
    if(dynamic_cast<QMainWindow* >(this)){
      settings.setValue("state", dynamic_cast<QMainWindow* >(this)->saveState());
    }
    settings.endGroup();
  }
  QString name_;
  QPoint pos_;
  QSize size_;
  bool visible_;
  bool default_;
};

template<class BASE, typename T1, typename T2>
class WidgetStateSaver<BASE,T1,T2>: public BASE
{
public:
  WidgetStateSaver(const QString& name, QWidget * parent ,T1 t1, T2 t2):
  BASE(parent,t1,t2),
  name_(name),
  pos_(),
  size_(),
  visible_(true),
  default_(false)
  {
    set_state_();
  }
  
  WidgetStateSaver(const QString& name,QPoint pos, QSize size, bool visible,
                   QWidget * parent ,T1 t1, T2 t2):
  BASE(parent,t1,t2),
  name_(name),
  pos_(pos),
  size_(size),
  visible_(visible),
  default_(true)
  {
    set_state_();
  }

  void ResetStateSaver(const QString& name)
  {
    name_=name;
    set_state_();
  }

protected:
  
  void set_state_()
  {
    QSettings settings;
    settings.beginGroup(name_);
    if(settings.contains("pos") && ! WidgetStateSaverConfig::SettingsIgnored()){
      BASE::move(settings.value("pos").toPoint());
      BASE::resize(settings.value("size").toSize());
      //BASE::setVisible(settings.value("visible").toBool());
      if(dynamic_cast<QMainWindow* >(this)){
        dynamic_cast<QMainWindow* >(this)->restoreState(settings.value("state").toByteArray());
      }
    }else if(default_){
      BASE::move(pos_);
      BASE::resize(size_);
      //BASE::setVisible(visible_);      
    }
    settings.endGroup();
  }
  
  virtual ~WidgetStateSaver()
  {
    QSettings settings;
    settings.beginGroup(name_);
    settings.setValue("pos", BASE::pos());
    settings.setValue("size", BASE::size());
    //settings.setValue("visible", /*BASE::isVisible()*/visible_);
    if(dynamic_cast<QMainWindow* >(this)){
      settings.setValue("state", dynamic_cast<QMainWindow* >(this)->saveState());
    }
    settings.endGroup();
  }
  QString name_;
  QPoint pos_;
  QSize size_;
  bool  visible_;
  bool default_;
};
template<class BASE, typename T1>
class WidgetStateSaver<BASE,T1>: public BASE
{
public:
  WidgetStateSaver(const QString& name, QWidget * parent ,T1 t1):
  BASE(parent,t1),
  name_(name),
  pos_(),
  size_(),
  visible_(true),
  default_(false)
  {
    set_state_();
  }
  
  WidgetStateSaver(const QString& name,QPoint pos, QSize size, bool visible,
                   QWidget * parent ,T1 t1):
  BASE(parent,t1),
  name_(name),
  pos_(pos),
  size_(size),
  visible_(visible),
  default_(true)
  {
    set_state_();
  }

  void ResetStateSaver(const QString& name)
  {
    name_=name;
    set_state_();
  }

protected:
  
  void set_state_()
  {
    QSettings settings;
    settings.beginGroup(name_);
    if(settings.contains("pos") && ! WidgetStateSaverConfig::SettingsIgnored()){
      BASE::move(settings.value("pos").toPoint());
      BASE::resize(settings.value("size").toSize());
      //BASE::setVisible(settings.value("visible").toBool());
      if(dynamic_cast<QMainWindow* >(this)){
        dynamic_cast<QMainWindow* >(this)->restoreState(settings.value("state").toByteArray());
      }
    }else if(default_){
      BASE::move(pos_);
      BASE::resize(size_);
      //BASE::setVisible(visible_);
    }
    settings.endGroup();
  }
  
  virtual ~WidgetStateSaver()
  {
    QSettings settings;
    settings.beginGroup(name_);
    settings.setValue("pos", BASE::pos());
    settings.setValue("size", BASE::size());
    //settings.setValue("visible", /*BASE::isVisible()*/visible_);
    if(dynamic_cast<QMainWindow* >(this)){
      settings.setValue("state", dynamic_cast<QMainWindow* >(this)->saveState());
    }
    settings.endGroup();
  }
  QString name_;
  QPoint pos_;
  QSize size_;
  bool   visible_;
  bool default_;
};
template<class BASE>
class WidgetStateSaver<BASE>: public BASE
{
public:
  WidgetStateSaver(const QString& name, QWidget * parent ):
  BASE(parent),
  name_(name),
  pos_(),
  size_(),
  visible_(true),
  default_(false)
  {
    set_state_();
  }
  
  WidgetStateSaver(const QString& name,QPoint pos, QSize size,
                   QWidget* parent):
  BASE(parent),
  name_(name),
  pos_(pos),
  size_(size),
  visible_(true),
  default_(true)
  {
    set_state_();
  }

  void ResetStateSaver(const QString& name)
  {
    name_=name;
    set_state_();
  }

protected:
  
  void set_state_()
  {
    QSettings settings;
    settings.beginGroup(name_);
    if(settings.contains("pos") && ! WidgetStateSaverConfig::SettingsIgnored()){
      BASE::move(settings.value("pos").toPoint());
      BASE::resize(settings.value("size").toSize());
      //BASE::setVisible(settings.value("visible").toBool());
      if(dynamic_cast<QMainWindow* >(this)){
        dynamic_cast<QMainWindow* >(this)->restoreState(settings.value("state").toByteArray());
      }
    }else if(default_){
      BASE::move(pos_);
      BASE::resize(size_);
      //BASE::setVisible(visible_);
    }
    settings.endGroup();
  }
  
  virtual ~WidgetStateSaver()
  {
    QSettings settings;
    settings.beginGroup(name_);
    settings.setValue("pos", BASE::pos());
    settings.setValue("size", BASE::size());
    //settings.setValue("visible", /*BASE::isVisible()*/visible_);
    if(dynamic_cast<QMainWindow* >(this)){
      settings.setValue("state", dynamic_cast<QMainWindow* >(this)->saveState());
    }
    settings.endGroup();
  }
  QString name_;
  QPoint pos_;
  QSize size_;
  bool visible_;
  bool default_;
};

}}//ns

#endif /*WIDGET_STATE_SAVER_HH_*/
