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
  DataViewer proxy

  Authors: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMG_GUI_DATA_VIEWER_PROXY_H
#define IMG_GUI_DATA_VIEWER_PROXY_H

#include "sip_handler.hh"

#include <boost/shared_ptr.hpp>
#include <QPointer>

#include <ost/img/point.hh>
#include <ost/img/extent.hh>
#include <ost/img/normalizer_impl.hh>

#include <ost/gui/data_viewer/overlay_base_fw.hh>
#include <ost/gui/data_viewer/overlay_manager_fw.hh>
#include <ost/gui/data_viewer/data_viewer.hh>
#include <ost/gui/module_config.hh>

class QWidget;
class QString;

namespace ost { namespace img { namespace gui {

// fw decl
class DataViewer;

class TEMPLATE_EXPORT DataViewerProxy  {
public:
  DataViewerProxy(DataViewer* v);

  // data viewer interface

  void SetData(const ImageHandle& data);
  void SetName(const String& name);

  NormalizerPtr GetNormalizer() const;
  void Renormalize();
  void Recenter();
  void UpdateView();

  Extent GetSelection() const;
  void SetSelection(const Extent& selection);
  void SetSlab(int slab);
  int GetSlab() const;
  void SetZoomScale(Real zoomscale);
  Real GetZoomScale() const;
  void SetViewerMin(Real min);
  Real GetViewerMin() const;
  void SetViewerMax(Real max);
  Real GetViewerMax() const;
  void SetGamma(Real gamma);
  Real GetGamma() const;
  void SetInvert(bool invert);
  bool GetInvert() const;
  void SetOffset(const geom::Vec2& offset);
  geom::Vec2 GetOffset() const;

  int AddOverlay(const OverlayPtr& ov, bool make_active=true);
  OverlayManagerPtr GetOverlayManager() const;
  void ClearOverlays();
  void ClearActiveOverlay();

  void AddDockWidget(QWidget* w, const QString& name, bool shown=true, int area=2);
  void RemoveDockWidget(QWidget* w);

  void SetAntialiasing(bool f);

  // ensure that data viewer is visible
  void Show();

  //void ShowClickedPosition(bool show=true);
  //Vec3 GetClickedPosition();

  void Hide();  
  object GetButtonGroup();
  object GetQObject();
protected:
  QPointer<DataViewer> viewer_ptr_;
  DataViewer* viewer() const;
};

typedef boost::shared_ptr<DataViewerProxy> DataViewerProxyPtr;

}}}  //ns

#endif


