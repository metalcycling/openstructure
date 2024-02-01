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
#ifndef OST_GUI_GL_CANVAS_HH
#define OST_GUI_GL_CANVAS_HH

#include <ost/gfx/glwin_base.hh>
#include <ost/gui/tools/tool.hh>
#include <ost/gfx/input.hh>

// Qt includes must come last
#include <QOpenGLWindow>
#include <QSurfaceFormat>

// forward declaration
class QResizeEvent;
class QOpenGLFramebufferObject;
class QOpenGLContext;
class QOffscreenSurface;

namespace ost { namespace gui {

class GLCanvas : public QOpenGLWindow, public gfx::GLWinBase {  
  Q_OBJECT;

public:
  GLCanvas();
  virtual ~GLCanvas();

  // gfx::GLWinBase interface
  virtual void MakeActive();
  virtual void DoRefresh() {this->update(); }
  virtual void StatusMessage(const String& m);
  virtual bool HasStereo() const {return format().stereo();};
  virtual bool HasMultisample() const {return format().samples() > 1;}

  // QSurfaceFormat for GL context setup
  void SetDefaultFormat();
  void SetStereoFormat();

  // central point for sending input to the gfx layer
  void OnTransform(gfx::InputCommand, int indx, 
                   gfx::TransformTarget, Real val);

  void SetTestMode(bool f);

  // Grab images from framebuffer and dump to disk
  virtual void Export(const String& fname, unsigned int width, 
                      unsigned int height, bool transparent);

  virtual void Export(const String& fname, unsigned int width, 
                      unsigned int height, int max_samples, bool transparent);

  virtual void Export(const String& fname, bool transparent);

signals:
  void CustomContextMenuRequested(const QPoint& point);

protected:

  // QOpenGLWindow interface
  virtual void initializeGL();
  virtual void paintGL();
  virtual void resizeGL(int w, int h);
  
  // event handling
  virtual void mouseMoveEvent(QMouseEvent* event);
  virtual void mousePressEvent(QMouseEvent* event);
  virtual void mouseReleaseEvent(QMouseEvent* event);
  virtual void mouseDoubleClickEvent(QMouseEvent* event);  

  virtual void keyPressEvent(QKeyEvent* event);
  virtual void keyReleaseEvent(QKeyEvent* event);
  virtual void wheelEvent(QWheelEvent* event);
  virtual bool event(QEvent* event);

private:

  bool IsToolEvent(QInputEvent* event) const;
  MouseEvent::Buttons TranslateButtons(Qt::MouseButtons buttons) const;
  void HandleMousePressEvent(QMouseEvent* event);
  void HandleMouseReleaseEvent(QMouseEvent* event);
  void HandleMouseMoveEvent(QMouseEvent* event);

  QPoint last_pos_;
  bool show_beacon_;  
  bool bench_flag_;

  // for image export
  bool offscreen_flag_;
  QOpenGLContext* offscreen_context_;
  QOffscreenSurface* offscreen_surface_;
  QOpenGLFramebufferObject* offscreen_fbo_;
};

}} // ns

#endif
