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
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMG_GUI_OVERLAY_BASE_H
#define IMG_GUI_OVERLAY_BASE_H


#include <ost/base.hh>

#include <ost/gui/module_config.hh>
#include "overlay_base_fw.hh"

#include <QObject>
#include <QEvent>
#include <QPainter>
#include <QMenu>
#include <QMouseEvent>
#include <QKeyEvent>
namespace ost { namespace img { namespace gui {

// fw decl
class DataViewerPanel;

// overlay base class
class DLLEXPORT_OST_GUI Overlay: public QObject
{
  Q_OBJECT;

public:
  //! initialization requires name
  Overlay(const String& name);

  //! set name
  void SetName (const String& name);
  
  //! retrieve name
  String GetName() const;

  //! call when redrawing is necessary
  /*!
    params: the drawing context, the parent data viewer, 
    as well as a flag to indicate wether this overlay is
    currently active
  */
  virtual void OnDraw(QPainter& pnt,  DataViewerPanel* dvp, bool is_active);

  //! mouse event handler
  /*!
    if an overlay is active, it receives the mouse events. if it handles
    the event for itself, it should return true, otherwise false should
    be returned to give the data viewer a chance to handle the event
    for itself
  */
  virtual bool OnMouseEvent(QMouseEvent* e,  DataViewerPanel* dvp, 
                            const QPoint& lastmouse);

  //! key event handler
  /*!
    see comments for the mouse events
  */
  virtual bool OnKeyEvent(QKeyEvent* e,  DataViewerPanel* dvp);

  //! return overlay specific pulldown menu
  /*
    this menu will be available in the viewer, under the name of
    the overlay. if this method returns zero (the default), then
    the viewer will ignore this menu
  */
  virtual QMenu* GetMenu();

  virtual void OnMenuEvent(QAction* e);

public slots:
  void OnMenuEventSlot(QAction* e);

signals:
  void SettingsChanged();
  void InfoTextChanged(const QString& t);

private:
  String name_;

  Overlay(const Overlay&) {}
  Overlay& operator=(const Overlay&) {return *this;}
};

}}}  //ns

#endif
