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
  Author: Andreas Schenk
*/

#include <ost/img/alg/fft.hh>   
#include "fft_panel.hh"


#include <QDebug>
#include <QInputDialog>


namespace ost { namespace img { namespace gui {
    
FFTPanel::FFTPanel(const ImageHandle& parent_data, QWidget* parent):
  DataViewerPanelBase(parent_data,parent),
  size_(std::min<int>(256,std::min<int>(parent_data.GetSize()[0],parent_data.GetSize()[1]))),
  parent_position_(parent_data.GetExtent().GetCenter()),
  fft_data_(CreateImage(Extent(Size(size_,size_),Point(0,0)),
                              COMPLEX,HALF_FREQUENCY)),
  parent_observer_(parent_data)
{
  QAction* a_set_size=new QAction("Set FFT size",this);
  popupmenu_->addAction(a_set_size);
  connect(a_set_size, SIGNAL(triggered()), this, SLOT(ShowSizeDialog()));
  DataViewerPanelBase::SetData(fft_data_);
  update_fft();
}

FFTPanel::~FFTPanel()
{
}


void FFTPanel::ObserverUpdate()
{
  UpdateView(true);
}

void FFTPanel::SetPosition(const Point& p)
{
  parent_position_=p;
  update_fft();
}

void FFTPanel::ShowSizeDialog()
{
  bool ok;
  #if QT_VERSION >= 0x040500
    int i = QInputDialog::getInt(this, "Set FFT size","FFT size", size_, 1, 
                                 std::min<int>(parent_observer_.GetObservedData().GetSize()[0],
                                               parent_observer_.GetObservedData().GetSize()[1]),
                                 1, &ok);
  #else
    int i = QInputDialog::getInteger(this, "Set FFT size","FFT size", size_, 1,
                                     std::min<int>(parent_observer_.GetObservedData().GetSize()[0],
                                                   parent_observer_.GetObservedData().GetSize()[1]),
                                     1, &ok);
  #endif
  if (ok){
    SetFFTSize(i);
  }
}

void FFTPanel::SetFFTSize(unsigned int size)
{
  size_=size;
  fft_data_.Reset(Extent(Size(size_,size_),Point(0,0)),
                  COMPLEX,HALF_FREQUENCY);
  update_fft();
}

unsigned int FFTPanel::GetFFTSize()
{
  return size_;
}

void FFTPanel::SetData(const ImageHandle& parent_data)
{
  parent_observer_=ParentDataObserver(parent_data);
}


void FFTPanel::update_fft()
{
  ImageHandle im=CreateImage(Extent(Size(size_,size_),parent_position_));
  im.Paste(parent_observer_.GetObservedData());
  im.ApplyIP(alg::FFT());
  fft_data_.Paste(im);
  UpdateView(true);
}


}}}  //ns

