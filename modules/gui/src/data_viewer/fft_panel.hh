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

#ifndef OST_FFT_PANEL_HH_
#define OST_FFT_PANEL_HH_
/*
  Author: Andreas Schenk
*/

#ifndef Q_MOC_RUN

#include <ost/gui/module_config.hh>
#include <ost/img/data_observer.hh>
#include "data_viewer_panel_base.hh"

#endif

namespace ost { namespace img { namespace gui {

class ParentDataObserver: public DataObserver
{
public:
  ParentDataObserver(const ImageHandle& d):DataObserver(d){}
  virtual void ObserverRelease()
  {
   // Noop as parent takes care of this.
  }
};

class DLLEXPORT_OST_GUI FFTPanel: public DataViewerPanelBase
{
Q_OBJECT

public:
  FFTPanel(const ImageHandle& parent_data, QWidget* parent);
  virtual ~FFTPanel();
  void SetFFTSize(unsigned int size);
  unsigned int GetFFTSize();
  virtual void ObserverUpdate();
  virtual void SetData(const ImageHandle& parent_data);

public slots:
  void SetPosition(const Point& p);
  void ShowSizeDialog();

protected:
  void update_fft();
  unsigned int size_;
  Point parent_position_;
  ImageHandle fft_data_;
  ParentDataObserver parent_observer_;
};

}}}  //ns

#endif /*FFT_PANEL_HH_*/
