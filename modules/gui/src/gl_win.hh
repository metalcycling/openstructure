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
#ifndef OST_GUI_GL_WIN_HH
#define OST_GUI_GL_WIN_HH


#include <ost/gfx/scene_observer.hh>
#include <ost/gui/module_config.hh>
#include <ost/gui/widget.hh>

/*
  Authors: Ansgar Philippsen, Marco Biasini
 */

namespace ost { namespace gui {


class GLCanvas;
class ToolBar;
class Tool;
class SceneMenu;

class DLLEXPORT_OST_GUI GLWin: public Widget, public gfx::SceneObserver
{
  Q_OBJECT;
public:
  GLWin(QWidget* p, bool try_stereo=false);
  ~GLWin();
  void SetTestMode(bool f);

public:
  virtual bool Restore(const QString& prefix);
  virtual bool Save(const QString& prefix);

public slots:
  void StatusMessage(const String& m);
  void ActiveToolChanged(Tool* t);

private slots:
  virtual void RequestContextMenu(const QPoint& pos);
  
private:
  GLCanvas*  gl_canvas_;
  QWidget*   gl_canvas_widget_;
  ToolBar*   toolbar_;
  SceneMenu* scene_menu_;
};

}} // ns

#endif

