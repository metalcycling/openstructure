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


#include <ost/base.hh>

#include "argand.hh"

#include <QPainter>
#include <QResizeEvent>
namespace ost { namespace img { namespace gui {

Argand::Argand(const ImageHandle& data, QWidget* p):
  QWidget(p),
  DataObserver(data),
  extent_(),
  current_(),
  has_data_(false),
  data_changed_(false),
  buffer_(new QPixmap(250,250)),
  sf_(1.0)
{
  setMinimumSize(100,100);
}

Argand::~Argand()
{
  delete buffer_;
}

void Argand::SetCurrentPixel(const Point& p)
{
 current_=p;
 update();  
}

void Argand::SetExtent(const Extent& e)
{
  has_data_=true;
  data_changed_=true;
  extent_=e;
  update();
}

void Argand::ClearExtent()
{
  has_data_=false;
  update();
}

void Argand::paintEvent(QPaintEvent* e)
{
  if(!IsDataValid()) return;
  QPainter pnt(this);
  if(has_data_) {
    if(data_changed_){
      data_changed_=false;
      update_buffer();
    }
    pnt.drawPixmap(QPoint(0,0),*buffer_);
  }

  if(extent_.Contains(current_) && GetObservedData().GetType()==COMPLEX){
    QPoint p0(size().width()/2,size().height()/2);
    Complex v = GetObservedData().GetComplex(current_);
    pnt.setPen(QPen(QColor::fromHsvF(fmod(std::arg(v)/(2*M_PI)+0.5,1.0),
                                     1,1.0),2));
    QPoint p((int)floor(v.real()*sf_),
             (int)floor(v.imag()*sf_));
    pnt.drawLine(p0,p0+p);
  }
}  

void Argand::update_buffer()
{
  if(!IsDataValid()) return;
  buffer_->fill();
  QPainter pnt(buffer_);
  QPoint p0(size().width()/2,size().height()/2);
  pnt.setPen(QPen(QColor(0,0,0)));
  pnt.drawLine(p0.x(),0,p0.x(),size().height());
  pnt.drawLine(0,p0.y(),size().width(),p0.y());
  if (!has_data_ || GetObservedData().GetType()!=COMPLEX) return;
  
  Real maxlen = 0.0;
  std::vector<Complex> data_list;
  for(ExtentIterator it(extent_); !it.AtEnd(); ++it) {
    Complex v = GetObservedData().GetComplex(it);
    pnt.setPen(QPen(QColor::fromHsvF(fmod(std::arg(v)/(2*M_PI)+0.5,1.0),1,0.8)));
    Real l = std::abs(v);
    maxlen = std::max(l,maxlen);
    data_list.push_back(v);
  }

  if(maxlen == 0.0) {
    return;    // return to avoid zero division
  }

  int dim=std::min(size().width(),size().height());
  sf_=0.5*(Real)(dim-10)/maxlen;

  for (std::vector<Complex>::iterator it=data_list.begin();
       it!=data_list.end();++it) {
    pnt.setPen(QPen(QColor::fromHsvF(fmod(std::arg(*it)/(2*M_PI)+0.5,1.0),
                    1,0.8)));
    QPoint p((int)floor(it->real()*sf_),
             (int)floor(it->imag()*sf_));
    pnt.drawLine(p0,p0+p);
  }
}
 
void Argand::resizeEvent(QResizeEvent* e)
{
  delete buffer_;
  buffer_=new QPixmap(e->size());
  data_changed_=true;
  update();
}

 
void Argand::ObserverUpdate()
{
  data_changed_=true;
  update();
}

void Argand::ObserverUpdate(const Extent& e)
{
  data_changed_=true;
  update();
}

void Argand::ObserverUpdate(const Point& p)
{
  data_changed_=true;
  update();
}

void Argand::ObserverRelease()
{
}

}}}  //ns
