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
  Authors: Ansgar Philippsen, Andreas Schenk, Jeff Lovelace
*/

#include <ost/gui/gosty_app.hh>
#include <cassert>
#include <iostream>
//#include <cstring>
//#include <sstream>


//#include <ost/io/io_manager.hh>
#include <ost/message.hh>
//#include <ost/img/alg/fft.hh>   
#include <ost/img/alg/norm.hh>   
#include <ost/gui/dock_widget.hh>

#include <ost/gui/perspective.hh>

#include "data_viewer.hh"
#include "data_viewer_panel.hh"
#include "overlay_base.hh"
#include "overlay_manager.hh"

#include <QLabel>
#include <QStatusBar>
#include <QMouseEvent>
#include <QVariant>
#include <QMenu>
#include <QMenuBar>
#include <QApplication>
#include <QDesktopWidget>
namespace ost { namespace img { namespace gui {

namespace {
int ipow(int base, unsigned int exponent){
  int result=1;
  for(unsigned int i=0;i<exponent;++i){
    result*=base;
  }
  return result;
}
}

DataViewer::DataViewer(QWidget* p, const ImageHandle& data, const QString& name):
  ost::gui::MainWindow(p),
  name_(name),
  panel_(new DataViewerPanel(data,this)),
  ov_manager_(new OverlayManager(this)),
  ov_manager_gui_(new OverlayManagerGUI(this, ov_manager_)),
  info_(new InfoPanel(this)),
  argand_(new Argand(data,this)),
  fft_(new FFTPanel(data,this)),
  lastmouse_()
{
  connect(ov_manager_gui_,SIGNAL(SettingsChanged()),this,SLOT(UpdateView()));
  setWindowTitle("OpenStructure Data Viewer");
  build(data);
  setFocusProxy(panel_);
}

DataViewer::~DataViewer()
{
  QMenuBar* mainMenu = ost::gui::GostyApp::Instance()->GetPerspective()->GetMenuBar();
  mainMenu->removeAction(WindowMenu()->menuAction());
}

void DataViewer::SetData(const ImageHandle& d)
{
  panel_->SetData(d);
  fft_->SetData(d);
  info_->SetImageInfo(d);
}

const ImageHandle& DataViewer::GetData() const{
  return panel_->GetObservedData();
}

NormalizerPtr DataViewer::GetNormalizer() const 
{
  assert(panel_);
  return (panel_->GetNormalizer());
}

//! re-apply normalization using current normalizer
void DataViewer::Renormalize() 
{
  assert(panel_);
  panel_->Renormalize();
}

void DataViewer::Recenter() 
{
  assert(panel_);
  panel_->Recenter();
}

//! update view
void DataViewer::UpdateView() 
{
  assert(panel_);
  panel_->UpdateView(false);
}

Extent DataViewer::GetSelection() const
{
  assert(panel_);
  return panel_->GetSelection();
}

void DataViewer::SetSelection(const Extent& selection)
{
  assert(panel_);
  panel_->SetSelection(selection);
}

void DataViewer::SetName(const String& name)
{
  setWindowTitle(QString::fromStdString(name));
}

int DataViewer::GetSlab() const 
{
  return panel_->GetSlab();
}

void DataViewer::SetZoomScale(Real zoomscale)
{
  panel_->SetZoomScale(zoomscale);
}

Real DataViewer::GetZoomScale() const
{
  return panel_->GetZoomScale();
}

void DataViewer::SetViewerMin(Real min)
{
  panel_->SetViewerMin(min);
}

Real DataViewer::GetViewerMin() const
{
  return panel_->GetViewerMin();
}

void DataViewer::SetViewerMax(Real max)
{
  panel_->SetViewerMax(max);
}

Real DataViewer::GetViewerMax() const
{
  return panel_->GetViewerMax();
}

void DataViewer::SetGamma(Real gamma)
{
  panel_->SetGamma(gamma);
}

Real DataViewer::GetGamma() const
{
  return panel_->GetGamma();
}

void DataViewer::SetInvert(bool invert)
{
  panel_->SetInvert(invert);
}

bool DataViewer::GetInvert() const
{
  return panel_->GetInvert();
}

void DataViewer::SetOffset(const geom::Vec2& offset)
{
  panel_->SetOffset(offset);
}

geom::Vec2 DataViewer::GetOffset() const
{
  return panel_->GetOffset();
}

QButtonGroup* DataViewer::GetButtonGroup(void)
{
  return ov_manager_gui_->GetButtonGroup();
}

int DataViewer::AddOverlay(const OverlayPtr& ov, bool make_active)
{
  int retval= ov_manager_->AddOverlay(ov,make_active);
  if(panel_) panel_->UpdateView(false);
  return retval;
}

void DataViewer::ClearOverlays()
{
  ov_manager_->DeleteAllOverlays();
}

void DataViewer::ClearActiveOverlay()
{
  ov_manager_->DeleteActiveOverlay();
}

OverlayManagerPtr DataViewer::GetOverlayManager() const
{
  return ov_manager_;
}

void DataViewer::SetSlab(int slab)
{
  panel_->SetSlab(slab);
}

void DataViewer::OnSlabChange(int slab)
{
  slablabel_->setText(QString::number(slab));
}

void DataViewer::OnZoomChange(int zoomlevel)
{
  if(zoomlevel>=0){
    zoomlabel_->setText(QString::number(ipow(2,zoomlevel))+":1");	
  }else{
    zoomlabel_->setText("1:"+QString::number(ipow(2,-zoomlevel)));	
  }
}


void DataViewer::OnPanelMouseEvent(QMouseEvent* event)
{
  geom::Vec2 mpos = panel_->WinToFracPoint(event->pos());
  Point mposi(mpos);
  bool lmb=event->buttons() == Qt::LeftButton;
  bool cmd=event->modifiers()&Qt::ControlModifier;

  if(panel_->GetObservedData().GetType()==REAL) {
    Real val=panel_->GetObservedData().GetReal(mposi);
    info_->SetMousePoint(mpos,val);
    if(lmb && !cmd) {
      info_->SetClickPoint(mpos,val);
    }
  } else {
    Complex val=panel_->GetObservedData().GetComplex(mposi);
    info_->SetMousePoint(mpos,val);
    if(lmb && !cmd) {
      info_->SetClickPoint(mpos,val);
    }
  }
}


void DataViewer::SetAntialiasing(bool f)
{
  panel_->SetAntialiasing(f);
}

//////////////////////////
// private methods

void DataViewer::build(const ImageHandle& data)
{
  build_menu();

  setAnimated(false);
  setCentralWidget(panel_);
  panel_->installEventFilter(this);
  show();
  connect(panel_,SIGNAL(zoomed(int)),SLOT(OnZoomChange(int)));
  connect(panel_,SIGNAL(slabChanged(int)),SLOT(OnSlabChange(int)));
  zoomlabel_=new QLabel("1:1");
  slablabel_=new QLabel("0");
  statusBar()->addWidget(zoomlabel_);
  statusBar()->addWidget(slablabel_);

  AddDockWidget(ov_manager_gui_,"Overlays",true,2);
  info_->SetImageInfo(data);
  AddDockWidget(info_,"Info",true,2);
  AddDockWidget(argand_,"Argand",false,1);
  AddDockWidget(fft_,"FFT",false,1);

  connect(panel_,SIGNAL(selected(const Extent&)),argand_,SLOT(SetExtent(const Extent&)));
  connect(panel_,SIGNAL(deselected()),argand_,SLOT(ClearExtent()));

  connect(panel_,SIGNAL(selected(const Extent&)),info_,SLOT(SetSelection(const Extent&)));
  connect(panel_,SIGNAL(deselected()),info_,SLOT(ClearSelection()));
  connect(panel_,SIGNAL(released()),this,SIGNAL(released()));
  
  if(!parentWidget()) {
    resize(QApplication::desktop()->availableGeometry().adjusted(20,20,-20,-20).size());
  } else {
    resize(1100,800);
  }
}

void DataViewer::build_menu()
{
QMenuBar* mainMenu = ost::gui::GostyApp::Instance()->GetPerspective()->GetMenuBar();
mainMenu->addMenu(WindowMenu());
}

 void DataViewer::show_current_window_menu(void)
   {
   // If the Panel menu for this window is already visible do nothing
   if (!WindowMenu()->isVisible())
     {
     // We need to hide all instances of panels
     // Get the Main Menu Bar pointer from OST
     QMenuBar* mainMenu = ost::gui::GostyApp::Instance()->GetPerspective()->GetMenuBar();
     // Get a list of all the main level menus
     QList<QAction*> ActionList=mainMenu->actions();
     // Any menus named panels are hidden
     for (int i=0;i<ActionList.count();i++)
       {
       if (ActionList[i]->text()=="Panels")
         {
         if (ActionList[i]->isVisible())
           {
           ActionList[i]->setVisible(false);
           }
         }
       }
     // Show the panel menu for this object
     WindowMenu()->menuAction()->setVisible(true);
     }
   }

bool DataViewer::eventFilter(QObject * object, QEvent *event)
{
  if (object == panel_) {
    switch (event->type()){
    case QEvent::KeyPress:
      {
        QKeyEvent *key_event = static_cast<QKeyEvent*>(event);
        if(key_event->modifiers() & Qt::ControlModifier) {
          ov_manager_->OnKeyEvent(key_event,panel_);
          return true;
        }
      }
      break;
    case QEvent::FocusIn:
      {
        // Show panels menu when data viewer gets focus
        show_current_window_menu();
        return true;
      }
      break;
    case QEvent::Show:
      {
        // Show panels menu when data viewer gets focus
        show_current_window_menu();
        return true;
      }
    case QEvent::MouseButtonDblClick:
    case QEvent::MouseButtonPress:
    case QEvent::MouseButtonRelease:
    case QEvent::MouseMove:
      {
        QMouseEvent *mouse_event = static_cast<QMouseEvent*>(event);
        if((mouse_event->modifiers() & Qt::ShiftModifier) && fft_->isVisible() )
	{
          fft_->SetPosition(panel_->WinToPoint(mouse_event->pos()));
        }
        OnPanelMouseEvent(mouse_event);
        if(mouse_event->modifiers() & Qt::ControlModifier) {
          ov_manager_->OnMouseEvent(mouse_event,panel_,lastmouse_);
          lastmouse_=mouse_event->pos();
          return true;
        }else{
          lastmouse_=mouse_event->pos();
        }
      }
      break;
    case QEvent::Paint:
      {
        QPaintEvent *paint_event = static_cast<QPaintEvent*>(event);
        panel_->paintEvent(paint_event);
        QPainter painter(panel_);
        ov_manager_->OnDraw(painter,panel_);
        return true;
      }
      break;
    default:
      break;
    }
    return false;
  }else{
    return MainWindow::eventFilter(object, event);
  }
}

}}}  //ns
