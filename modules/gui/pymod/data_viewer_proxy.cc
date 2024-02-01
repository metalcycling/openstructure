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


#include "data_viewer_proxy.hh"
#include <ost/message.hh>

namespace ost { namespace img { namespace gui {

DataViewerProxy::DataViewerProxy(DataViewer* v):
  viewer_ptr_(v)
{
}


NormalizerPtr DataViewerProxy::GetNormalizer() const
{
  return viewer()->GetNormalizer();
}

void DataViewerProxy::Renormalize()
{
  viewer()->Renormalize();
}

void DataViewerProxy::UpdateView()
{
  viewer()->UpdateView();
}

void DataViewerProxy::Recenter()
{
  viewer()->Recenter();
}

Extent DataViewerProxy::GetSelection() const
{
  return viewer()->GetSelection();
}
void DataViewerProxy::SetSelection(const Extent& selection)
{
  viewer()->SetSelection(selection);
}

void DataViewerProxy::SetSlab(int slab)
{
  viewer()->SetSlab(slab);
}

int DataViewerProxy::GetSlab() const
{
  return viewer()->GetSlab();
}

void DataViewerProxy::SetZoomScale(Real zoomscale)
{
  viewer()->SetZoomScale(zoomscale);
}

Real DataViewerProxy::GetZoomScale() const
{
  return viewer()->GetZoomScale();
}

void DataViewerProxy::SetViewerMin(Real min)
{
  viewer()->SetViewerMin(min);
}

Real DataViewerProxy::GetViewerMin() const
{
  return viewer()->GetViewerMin();
}

void DataViewerProxy::SetViewerMax(Real max)
{
  viewer()->SetViewerMax(max);
}

Real DataViewerProxy::GetViewerMax() const
{
  return viewer()->GetViewerMax();
}

void DataViewerProxy::SetGamma(Real gamma)
{
  viewer()->SetGamma(gamma);
}

Real DataViewerProxy::GetGamma() const
{
  return viewer()->GetGamma();
}

void DataViewerProxy::SetInvert(bool invert)
{
  viewer()->SetInvert(invert);
}

bool DataViewerProxy::GetInvert() const
{
  return viewer()->GetInvert();
}

void DataViewerProxy::SetOffset(const geom::Vec2& offset)
{
  viewer()->SetOffset(offset);
}

geom::Vec2 DataViewerProxy::GetOffset() const
{
  return viewer()->GetOffset();
}


void DataViewerProxy::SetData(const ImageHandle& d)
{
  viewer()->SetData(d);
}

void DataViewerProxy::SetName(const String& name)
{
  viewer()->SetName(name);
}

int DataViewerProxy::AddOverlay(const OverlayPtr& ov, bool make_active)
{
  return viewer()->AddOverlay(ov,make_active);
}

void DataViewerProxy::ClearOverlays()
{
  viewer()->ClearOverlays();
}

void DataViewerProxy::ClearActiveOverlay()
{
  viewer()->ClearActiveOverlay();
}

OverlayManagerPtr DataViewerProxy::GetOverlayManager() const
{
  return viewer()->GetOverlayManager();
}


void DataViewerProxy::AddDockWidget(QWidget* w, const QString& name, bool shown, int area)
{
  viewer()->AddDockWidget(w,name, shown,area);
}

void DataViewerProxy::RemoveDockWidget(QWidget* w)
{
  viewer()->RemoveDockWidget(w);
}

void DataViewerProxy::SetAntialiasing(bool f)
{
  viewer()->SetAntialiasing(f);
}

void DataViewerProxy::Show()
{
  viewer()->show();
}

void DataViewerProxy::Hide() 
{
  viewer()->hide();
}

object DataViewerProxy::GetButtonGroup()
{
  return ost::gui::get_py_qobject<QButtonGroup>(viewer()->GetButtonGroup());
}

object DataViewerProxy::GetQObject(){
  return ost::gui::get_py_qobject<DataViewer>(viewer());
}

DataViewer* DataViewerProxy::viewer() const
{
  if(!viewer_ptr_){
     throw Error("DataViewer no longer exists");
  }
  return viewer_ptr_;
}
}}}  //ns
