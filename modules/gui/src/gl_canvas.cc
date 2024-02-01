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

#include <chrono>

#include "gl_canvas.hh"

#include <ost/log.hh>
#include <ost/gfx/scene.hh>
#include <ost/gui/tools/tool_manager.hh>
#include <ost/gui/perspective.hh>
#include <ost/gui/gosty_app.hh>

#include <QResizeEvent>
#include <QMouseEvent>
#include <QImage>
#include <QString>
#include <QOpenGLFramebufferObject>
#include <QOffscreenSurface>
#include <QOpenGLContext>
#include <QOpenGLFunctions>


namespace ost { namespace gui {

ost::gui::GLCanvas::GLCanvas(): QOpenGLWindow(),
                                last_pos_(),
                                show_beacon_(false),
                                bench_flag_(false),
                                offscreen_flag_(false),
                                offscreen_context_(NULL),
                                offscreen_surface_(NULL),
                                offscreen_fbo_(NULL) { 
  LOG_DEBUG("GLCanvas::registering with scene");
  gfx::Scene::Instance().Register(this);  
}

ost::gui::GLCanvas::~GLCanvas() {
  gfx::Scene::Instance().Unregister(this);
  if(offscreen_fbo_ != NULL) {
    // all other offscreen rendering objects are also != NULL
    // cleanup done as in QT threadrenderer example, not sure whether
    // offscreen context must be made current to delete offscreen buffer...
    offscreen_context_->makeCurrent(offscreen_surface_);
    delete offscreen_fbo_;
    offscreen_context_->doneCurrent();
    delete offscreen_context_;
    delete offscreen_surface_;
    offscreen_flag_ = false;
  }
}

void GLCanvas::MakeActive() {
  if(offscreen_flag_) {
    offscreen_context_->makeCurrent(offscreen_surface_);
  } else {
    this->makeCurrent();
  }
}

void GLCanvas::StatusMessage(const String& m) {
  // This Window can also be displayed without a full blown GostyApp.
  // We therefore only feed the message into Gosty if there's already
  // a valid instance.
  if(GostyApp::ValidInstance()) {
    GostyApp::Instance()->GetPerspective()->StatusMessage(m);
  }
}

void GLCanvas::SetDefaultFormat() {
  QSurfaceFormat f = QSurfaceFormat::defaultFormat();
  f.setRedBufferSize(8);
  f.setGreenBufferSize(8);
  f.setBlueBufferSize(8);
  f.setAlphaBufferSize(8);
  f.setDepthBufferSize(24);
  this->setFormat(f);
}

void GLCanvas::SetStereoFormat() {
  QSurfaceFormat f = QSurfaceFormat::defaultFormat();
  f.setRedBufferSize(8);
  f.setGreenBufferSize(8);
  f.setBlueBufferSize(8);
  // QOpenGLWindow seems to dislike alphabuffer in stereo rendering...
  //f.setAlphaBufferSize(8);
  f.setDepthBufferSize(24);
  f.setStereo(true);
  this->setFormat(f);
}

void GLCanvas::initializeGL() {
  LOG_DEBUG("GLCanvas::initializeGL()");
  gfx::Scene::Instance().InitGL();    
}

void GLCanvas::paintGL() {

  // static variables for benchmarking
  static int benchmark_count = 0;
  static std::chrono::steady_clock::time_point timer;

  if(bench_flag_) {
    ++benchmark_count;
    if(benchmark_count == 20) {
      benchmark_count = 0;
      int ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - timer).count();
      Real fps = Real(20) / (Real(ms) * 0.001);
      LOG_ERROR(fps << " fps");
      timer = std::chrono::steady_clock::now();
    }
    this->update();
  } 

  gfx::Scene::Instance().RenderGL(); 
  Tool* tool=ToolManager::Instance().GetActiveTool();
  if (tool) {
    tool->RenderGL();
  }
}

void GLCanvas::resizeGL(int w, int h) {
  gfx::Scene::Instance().Resize(w, h);
}

void GLCanvas::mouseMoveEvent(QMouseEvent* event) {
  
  if (!(show_beacon_ || event->buttons())) {
    return;
  }

  if (this->IsToolEvent(event)) {
    if (ToolManager::Instance().GetActiveTool()) {
      MouseEvent mouse_event(this->TranslateButtons(event->buttons()), 
                             last_pos_, QPoint(event->x(), event->y()),
                             event->modifiers());
      ToolManager::Instance().GetActiveTool()->MouseMove(mouse_event);
      event->accept();      
    }
  } else {
    this->HandleMouseMoveEvent(event);
  }
  last_pos_=QPoint(event->x(), event->y());
  this->DoRefresh();
}

void GLCanvas::mousePressEvent(QMouseEvent* event) {
  
  if (this->IsToolEvent(event)) {
    if (ToolManager::Instance().GetActiveTool()) {
      MouseEvent mouse_event(this->TranslateButtons(event->button()), 
                             QPoint(event->x(), event->y()),
                             QPoint(event->x(), event->y()),
                             event->modifiers());
      ToolManager::Instance().GetActiveTool()->MousePress(mouse_event);
    }
  } else {
    this->HandleMousePressEvent(event);
  }
  last_pos_=QPoint(event->x(), event->y());
  this->DoRefresh();
}

void GLCanvas::mouseReleaseEvent(QMouseEvent* event) {
  
  if (this->IsToolEvent(event)) {
    if (ToolManager::Instance().GetActiveTool()) {
      MouseEvent mouse_event(this->TranslateButtons(event->button()), 
                             QPoint(event->x(), event->y()),
                             QPoint(event->x(), event->y()),
                             event->modifiers());
      ToolManager::Instance().GetActiveTool()->Click(mouse_event);
      event->accept();      
    } else {
      this->HandleMouseReleaseEvent(event);
    }
  }
  last_pos_=QPoint(event->x(), event->y());
  this->DoRefresh();
}

void GLCanvas::mouseDoubleClickEvent(QMouseEvent* event) {
  
  if (this->IsToolEvent(event)) {
    if (ToolManager::Instance().GetActiveTool()) {
      MouseEvent mouse_event(this->TranslateButtons(event->button()), 
                             QPoint(event->x(), event->y()),
                             QPoint(event->x(), event->y()),
                             event->modifiers());
      ToolManager::Instance().GetActiveTool()->DoubleClick(mouse_event);
      event->accept();
    }
  } else {
    gfx::Scene& scene=gfx::Scene::Instance();
    mol::AtomHandle atom=scene.PickAtom(event->x(), 
                                        size().height()-event->y()).second;
    if (atom.IsValid()) {
      scene.SetCenter(atom.GetPos());
    }
  }
  if (show_beacon_) {
    this->DoRefresh();
  }
}

void GLCanvas::OnTransform(gfx::InputCommand com, int indx, 
                           gfx::TransformTarget trg, Real val) {
  // does not request a redraw on purpose, in order to chain OnTransform calls
  ost::gfx::Scene::Instance().Apply(gfx::InputEvent(gfx::INPUT_DEVICE_MOUSE, com, indx,
                                                    trg, val*0.5),false);
}

void GLCanvas::SetTestMode(bool f) {
  gfx::Scene::Instance().SetTestMode(f);
}

void GLCanvas::Export(const String& fname, unsigned int width, 
                      unsigned int height, bool transparent) {
  this->Export(fname, width, height, 0, transparent);
}

void GLCanvas::Export(const String& fname, unsigned int width, 
                      unsigned int height, int max_samples, bool transparent) {

  // setup of context, surface, fbo etc are implemented as in the QT
  // threadrenderer example
  gfx::Viewport old_vp = gfx::Scene::Instance().GetViewport();

  if(old_vp.width == static_cast<int>(width) && 
     old_vp.height == static_cast<int>(height) && max_samples <= 0) {
    // just grab the framebuffer, no need for fancy offscreen rendering...
    this->Export(fname, transparent);
    return;
  }

  offscreen_flag_ = true;

  if(offscreen_context_ == NULL || offscreen_surface_ == NULL) {
    // Some GL implementations require that the currently bound context is
    // made non-current before we set up sharing, so we doneCurrent here
    // and makeCurrent down below while setting up our own context.
    this->context()->doneCurrent();
    offscreen_context_ = new QOpenGLContext();
    QSurfaceFormat f = this->context()->format();
    if(max_samples > 0) {
      f.setSamples(max_samples);
    } 
    offscreen_context_->setFormat(f);
    offscreen_context_->setShareContext(this->context());
    offscreen_context_->create();

    offscreen_surface_ = new QOffscreenSurface();
    offscreen_surface_->setFormat(offscreen_context_->format());
    offscreen_surface_->create();

    offscreen_context_->makeCurrent(offscreen_surface_);
    gfx::Scene::Instance().ContextSwitch();
    gfx::Scene::Instance().InitGL(false);
  } else {
    this->context()->doneCurrent();
    // The following line destroys and recreates the existing OpenGL context.
    // This is an ugly fix to avoid weird artefacts when calling the Export 
    // function the second time. On a low-end Laptop this requires ~5-10ms and
    // is therefore much faster than the rate-limiting saving of the image to
    // disk at the end of the function. Feel free to investigate and propose a
    // more efficient solution.
    offscreen_context_->create();
    offscreen_context_->makeCurrent(offscreen_surface_);
    gfx::Scene::Instance().ContextSwitch();
    // the following InitGL sets potentially changed glClearcolor etc
    // could be made more efficient...
    gfx::Scene::Instance().InitGL(false);
  }

  if(offscreen_fbo_ == NULL || 
     offscreen_fbo_->width() != static_cast<int>(width) || 
     offscreen_fbo_->height() != static_cast<int>(height)) {
    if(offscreen_fbo_ != NULL) {
      delete offscreen_fbo_;
    }
    QOpenGLFramebufferObjectFormat fbo_format;
    // the following flag is required for OpenGL depth testing
    fbo_format.setAttachment(QOpenGLFramebufferObject::CombinedDepthStencil);
    offscreen_fbo_ = new QOpenGLFramebufferObject(width, height, fbo_format);
  }

  offscreen_fbo_->bind();
  gfx::Scene::Instance().Resize(width, height);
  this->paintGL();

  if(!transparent) {
    gfx::Scene::Instance().SetAlphaBias(1.0);
  }

  offscreen_context_->functions()->glFlush();
  QImage image = offscreen_fbo_->toImage();
  offscreen_fbo_->release();

  if(!transparent) {
    gfx::Scene::Instance().SetAlphaBias(0.0);
  }

  image.save(QString(fname.c_str()));
  offscreen_flag_ = false;
  this->MakeActive();
  gfx::Scene::Instance().Resize(old_vp.width, old_vp.height);
  gfx::Scene::Instance().ContextSwitch();
}

void GLCanvas::Export(const String& fname, bool transparent) {

  if(!transparent) {
    gfx::Scene::Instance().SetAlphaBias(1.0);
  }
  QImage image = this->grabFramebuffer();
  if(!transparent) {
    gfx::Scene::Instance().SetAlphaBias(0.0);
  }
  image.save(QString(fname.c_str()));
}

bool GLCanvas::IsToolEvent(QInputEvent* event) const {
  return event->modifiers() & Qt::ControlModifier;
}

MouseEvent::Buttons GLCanvas::TranslateButtons(Qt::MouseButtons buttons) const {
  return MouseEvent::Buttons(buttons);
}

void GLCanvas::HandleMousePressEvent(QMouseEvent* event) {
  gfx::Scene& scene=gfx::Scene::Instance();
  scene.Pick(event->x(), scene.GetViewport().height-event->y(), 0);
  event->accept();

  // if its a right click, we emit the ContextMenu signal,
  // maybe someone wants to implement something interactive on top...
  if(event->buttons() == Qt::RightButton) {
    emit CustomContextMenuRequested(event->pos());
  } 
}

void GLCanvas::HandleMouseReleaseEvent(QMouseEvent* event) {
  event->accept();  
}

void GLCanvas::HandleMouseMoveEvent(QMouseEvent* event) {
  int indx=0;
  gfx::TransformTarget trg=gfx::TRANSFORM_VIEW;  

  if(show_beacon_) {
    ost::gfx::Scene::Instance().SetBeacon(event->x(),size().height()-event->y());
  }

  QPoint delta=QPoint(event->x(), event->y())-last_pos_;
  if (event->buttons() & Qt::LeftButton) {
    if (event->buttons() & Qt::MiddleButton) {
      if (event->modifiers() & Qt::ShiftModifier) {
        OnTransform(gfx::INPUT_COMMAND_SLABN,indx,
                    trg,0.5*static_cast<Real>(-delta.x()));
        OnTransform(gfx::INPUT_COMMAND_SLABF,indx,
                    trg,0.5*static_cast<Real>(delta.x()));

      } else {
        OnTransform(gfx::INPUT_COMMAND_ROTZ, indx, trg, 
                    static_cast<Real>(delta.x()+delta.y()));
      }

    } else {
      if (event->modifiers() & Qt::ShiftModifier) {
        OnTransform(gfx::INPUT_COMMAND_TRANSX,indx,
                    trg,static_cast<Real>(delta.x())); // translate along x
        OnTransform(gfx::INPUT_COMMAND_TRANSY,indx,
                    trg,static_cast<Real>(-delta.y())); // translate along y
      } else {
        OnTransform(gfx::INPUT_COMMAND_ROTY,indx,trg,
                    static_cast<Real>(delta.x())); // rotation around y
        OnTransform(gfx::INPUT_COMMAND_ROTX,indx,trg,
                    static_cast<Real>(delta.y())); // rotation around x
      }

    }
  } else if (event->buttons() & Qt::MiddleButton) {
      if (event->modifiers() & Qt::ShiftModifier) {
        OnTransform(gfx::INPUT_COMMAND_SLABN,indx,trg,
                    0.5*static_cast<Real>(-delta.x()));
        OnTransform(gfx::INPUT_COMMAND_SLABF,indx,trg,
                    0.5*static_cast<Real>(-delta.x()));
      } else {
        OnTransform(gfx::INPUT_COMMAND_TRANSZ,indx,
                    trg,static_cast<Real>(delta.y())); // translate along z
      }
  }
  event->accept();
}

void GLCanvas::keyPressEvent(QKeyEvent* event) {
  if(event->key()==Qt::Key_Space) {
    show_beacon_=true;
    gfx::Scene::Instance().SetBeacon(last_pos_.x(),size().height()-last_pos_.y());
    DoRefresh();
    setCursor(Qt::BlankCursor);
  }
  if((event->modifiers() & Qt::ControlModifier)) {
    // Ctrl pressed
    if(event->key()==Qt::Key_A) {
      OnTransform(gfx::INPUT_COMMAND_AUTOSLAB,0,gfx::TRANSFORM_VIEW,0.0);
      return;
    } else if(event->key()==Qt::Key_F) {
      OnTransform(gfx::INPUT_COMMAND_TOGGLE_FOG,0,gfx::TRANSFORM_VIEW,0.0);
      return;
    } else if(event->key()==Qt::Key_B) { 
      bench_flag_=!bench_flag_;
      return;
    } else if(event->key()==Qt::Key_1) {
      gfx::Scene::Instance().SetShadingMode("fallback");
      DoRefresh();
      return;
    } else if(event->key()==Qt::Key_2) {
      gfx::Scene::Instance().SetShadingMode("basic");
      DoRefresh();
      return;
    } else if(event->key()==Qt::Key_3) {
      gfx::Scene::Instance().SetShadingMode("default");
      DoRefresh();
      return;
    } else if(event->key()==Qt::Key_4) {
      gfx::Scene::Instance().SetShadingMode("hf");
      DoRefresh();
      return;
    } else if(event->key()==Qt::Key_5) {
      gfx::Scene::Instance().SetShadingMode("toon1");
      DoRefresh();
      return;
    } else if(event->key()==Qt::Key_6) {
      gfx::Scene::Instance().SetShadingMode("toon2");
      DoRefresh();
      return;
    } else if(event->key()==Qt::Key_0) {
      gfx::Scene::Instance().SetShadow(!gfx::Scene::Instance().GetShadow());
      DoRefresh();
      return;
    } else if(event->key()==Qt::Key_9) {
      gfx::Scene::Instance().SetAmbientOcclusion(!gfx::Scene::Instance().GetAmbientOcclusion());
      DoRefresh();
      return;
    } else if(event->key()==Qt::Key_Equal) {
      if(gfx::Scene::Instance().GetStereoMode()>0) {
        gfx::Scene::Instance().SetStereoMode(0);
      } else {
        gfx::Scene::Instance().SetStereoMode(1);
      }
      DoRefresh();
      return;
    }    
  }
  event->ignore();
}

void GLCanvas::keyReleaseEvent(QKeyEvent* event) {
  if(event->key()==Qt::Key_Space) {
    show_beacon_=false;
    gfx::Scene::Instance().SetBeaconOff();
    DoRefresh();
    setCursor(Qt::ArrowCursor);
    return;
  }
}

void GLCanvas::wheelEvent(QWheelEvent* event) {
  OnTransform(gfx::INPUT_COMMAND_TRANSZ,0,gfx::TRANSFORM_VIEW,
              0.2*static_cast<Real>(-event->angleDelta().y()));
}

bool GLCanvas::event(QEvent* event) {
  return QOpenGLWindow::event(event);
}

}} // ns
